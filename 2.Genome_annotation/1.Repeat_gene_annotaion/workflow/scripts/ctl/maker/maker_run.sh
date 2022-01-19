threads=$1
name=$2
module load maker/2.31.11
mpiexec -n $thread maker -base $name maker_bopts.ctl maker_exe.ctl maker_opts.ctl 1>maker.log 2>maker.err
