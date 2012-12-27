defmode='install'
arg=${1:-$defmode}

if [ ! -e logs ]
then
    mkdir logs/
fi

sudo /data/pyramid/bin/python setup.py build > logs/build.log
sudo /data/pyramid/bin/python setup.py $arg > logs/deploy.log
sudo chmod -R 775 *
sudo chown -R www-data *
sudo chgrp -R development *

if [ -e "logs/ptmworker.pid" ]
then
    kill `cat "logs/ptmworker.pid"`
    while [ -f "logs/ptmworker.pid" ]
    do
        sleep 1
    done
fi

if [ "$arg" == "develop" ]
then
    /data/pyramid/bin/pceleryd development.ini --logfile=logs/ptmworker.log --loglevel=DEBUG --pidfile=logs/ptmworker.pid > logs/ptmworker.start &
else
    /data/pyramid/bin/pceleryd production.ini --logfile=logs/ptmworker.log --loglevel=INFO --pidfile=logs/ptmworker.pid > logs/ptmworker.start &
fi
sudo apache2ctl restart

