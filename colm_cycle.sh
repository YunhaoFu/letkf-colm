#!/bin/sh
# Author: Yongjun Zheng, Sunday Apr 2, 2023

if [ $# -lt 4 ]; then
  echo "Usage: colm_cycle.sh start_YMDH:end_YMDH cycle_hours n_members forcing_dir da_executable obs_dir"
  exit
fi

start_YMDH=`echo $1 | awk '{split($1, f, ":"); print f[1]}'`
end_YMDH=`echo $1 | awk '{split($1, f, ":"); print f[2]}'`
[ -z "$end_YMDH" ] && echo "You must specify end_YMDH" && exit 1

cycle_hours=$2
n_members=$3
forcing_dir=$4
da_executable=${5:-}
obs_dir=${6:-} # please give observation path if you would like to do data assimiliation

if [ -z "$da_executable" ]; then
  with_da=0
elif [ -f $da_executable ]; then
  with_da=1
  [ -z "$obs_dir" ] && echo "You must specify obs_dir to do DA" && exit 1
else
  with_da=0
fi

########################################################################
# check restart files of members
########################################################################

########################################################################
# cycle loop
########################################################################
nx=720
ny=360
dt=900
nt=$((3600*cycle_hours/dt))
YMDH=${start_YMDH}
while [ ${YMDH} -le ${end_YMDH} ]
do
  # current cycle start time
  start_yr=${YMDH:0:4}
  start_jday=`date +%j -u -d "${YMDH:0:8}"`
  HH=${YMDH:8:2}
  start_sec=$((HH*3600))
  
  # current cycle end time
  YMDH=`date +%Y%m%d%H -u -d "${YMDH:0:8} ${YMDH:8:2} +${cycle_hours}hours"`
  end_yr=${YMDH:0:4}
  end_jday=`date +%j -u -d "${YMDH:0:8}"`
  HH=${YMDH:8:2}
  end_sec=$((HH*3600))
  if [ -d ${obs_dir} ]; then
    obs=${obs_dir}/${end_yr}/${end_yr}-${end_jday}-${end_sec}.nc
    nt_current=${nt}
    until [ -f ${obs} ]; do
      nt_current=$((nt_current+nt))

      YMDH=`date +%Y%m%d%H -u -d "${YMDH:0:8} ${YMDH:8:2} +${cycle_hours}hours"`
      end_yr=${YMDH:0:4}
      end_jday=`date +%j -u -d "${YMDH:0:8}"`
      HH=${YMDH:8:2}
      end_sec=$((HH*3600))
      
      obs=${obs_dir}/${end_yr}/${end_yr}-${end_jday}-${end_sec}.nc
    done
  fi

  # top path
  tdir=`pwd`
  
  ######################################################################
  # run each member for the period of a cycle
  ######################################################################
  echo "running CoLM from ${start_yr}-${start_jday}-${start_sec} to ${end_yr}-${end_jday}-${end_sec} in ${tdir}"
  for member in `seq -f '%03g' 1 ${n_members}`
  do
    # member path
    bkg_dir=${tdir}/bkg/${member}
    ana_dir=${tdir}/ana/${member}
    # temporal path under which a colm member will run
    tmp_dir=${tdir}/tmp/${member}
    if [ -d ${tmp_dir} ]; then
      rm -rf ${tmp_dir}/*
    else
      mkdir -p ${tmp_dir}
    fi
    cd ${tmp_dir}
    
    # generate job script
cat > ${tmp_dir}/colm_run.cmd << CMD
#!/bin/sh
#SBATCH --job-name=CoLM
#SBATCH --partition=Regular
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1

# Create an input parameter namelist file
cat > ${tmp_dir}/colm.nml << EOF;
&colm_nml
greenwich      =  .true.,
start_yr       =  ${start_yr},
start_jday     =  ${start_jday},
start_sec      =  ${start_sec},
fsurdat        = '',
flaidat        = '',
forcdir        = '${forcing_dir}',
frestart       = '${tmp_dir}/restart',
lon_points     =  ${nx},
lat_points     =  ${ny},
edgen          = 90.,
edgee          = 180.,
edges          = -90.,
edgew          = -180.,
deltim         =  ${dt},
mstep          =  ${nt_current},
ostep          =  ${nt},
/
EOF

with_da=${with_da}
if [ \${with_da} -eq 0 ]; then
  ln -sf ${bkg_dir}/restart-${start_yr}-${start_jday}-${start_sec}.nc ${tmp_dir}
else
  ln -sf ${ana_dir}/restart-${start_yr}-${start_jday}-${start_sec}.nc ${tmp_dir}
fi

NPROC=\${SLURM_NTASKS}
mpirun --mca pml ucx -mca osc ucx -n \${NPROC} ${tdir}/colm.exe
retval=\$?

if [ \${retval} -eq 0 ]; then
  mv ${tmp_dir}/restart-${end_yr}-${end_jday}-${end_sec}.nc ${bkg_dir}
  mv ${tmp_dir}/fldxy_instant-*.nc ${bkg_dir}
  mv ${tmp_dir}/fldxy_average-*.nc ${bkg_dir}

  touch ${tmp_dir}/.colm_done
fi
CMD
    
    # run
    sbatch ${tmp_dir}/colm_run.cmd
    sleep 2 # wait for two seconds before submitting another job to prevent from crashing the job system
  done
  
  ######################################################################
  # wait for the finishing of all members of a cycle
  ######################################################################
  n=0
  while [ ${n} -ne ${n_members} ]
  do
    sleep 30
    
    n=0
    for member in `seq -f '%03g' 1 ${n_members}`
    do
      # temporal path under which a colm member is running
      tmp_dir=${tdir}/tmp/${member}
      
      [ -f ${tmp_dir}/.colm_done ] && n=$((n+1))
    done
  done
  
  ######################################################################
  # do data assimilation
  ######################################################################
  if [ ${with_da} -eq 1 ]  & [ -f ${obs} ]; then
    n=0
    for member in `seq -f '%03g' 1 ${n_members}`
    do
      # member path
      bkg_dir=${tdir}/bkg/${member}
      ana_dir=${tdir}/ana/${member}
      # temporal path under which the data assimilation will run
      tmp_dir=${tdir}/tmp
      if [ -d ${tmp_dir} ]; then
        rm -rf ${tmp_dir}/*
      else
        mkdir -p ${tmp_dir}
      fi
      cd ${tmp_dir}
      
      bkg=${bkg_dir}/restart-${end_yr}-${end_jday}-${end_sec}.nc
      ana=${ana_dir}/restart-${end_yr}-${end_jday}-${end_sec}.nc
      if [ -f ${bkg} ]; then
        n=$((n+1))
        cp -f ${bkg} ${ana}
        ln -sf ${bkg} ${tmp_dir}/bkg.${member}
        ln -sf ${ana} ${tmp_dir}/ana.${member}
    done
    
    # do data assimilation if all the backgrounds are presented
    if [ ${n} -eq ${n_members} ]; then
      echo "doing DA at ${end_yr}-${end_jday}-${end_sec} in ${tdir}"
      ln -sf ${obs_dir}/${obs}.nc ${tmp_dir}/obs.nc
      
      # generate DA namelist
cat > ${tmp_dir}/input.nml << EOF
&description
idate = ${end_yr}, ${end_jday}, ${end_sec}
patch2ij_file = "./global_0.5_noreverse_invar.nc"
lon_points = ${nx} 
lat_points = ${ny} 
edges = -90.00
edgen =  90.00
edgew = -180.00
edgee =  180.00
ens_size = ${n_members}
max_y_delta = 5.0 
radius = 50000.0,5000.0
infl = 1.0 
save_diag = .true.
diag_dir = "./diag"
/
EOF

      # run DA
      ${da_executable}

      # wait for the finishing of DA
    else
      echo "Something wrong with the previous cycle, please check and fix, then restart at at ${start_yr}-${start_jday}-${start_sec}"
      exit
    fi
  fi
done
