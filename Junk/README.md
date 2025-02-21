# Parallel-Scientific-Computing

Repo for PSC courseworks

Use this for creating initial conditions
- python3 create_initial_conditions.py --final-time 60 --snapshots 1 --executable-name  ./step-0 --min-mass 0.003 --max-mass 0.01 --dt 1 --N 2


python3 create_initial_conditions.py --final-time 60 --snapshots 1 --executable-name  ./step-0 --min-mass 1 --max-mass 3 --dt 1 --N 10^C

dt = 0.1
n = 10

^^^^ reasonable results

then followed by chmod u+x ./step-1.sh to allow executable permissions

When creating a python config, manually change the executable to step-0, step-1 etc to check as parameters are different.


SUBMISSION 2 python script:
python3 create_initial_conditions.py --final-time 300 --snapshots 20 --executable-name  ./step-0 --min-
mass 0.003 --max-mass 0.01 --dt 1 --N 10







