defmode='install'
arg=${1:-$defmode}

if [ ! -e logs ]
then
    mkdir logs/
fi

echo "Deploying codebase..."
sudo /data/pyramid/bin/python setup.py build > logs/build.log
sudo /data/pyramid/bin/python setup.py $arg > logs/deploy.log
sudo chmod -R 775 *
sudo chown -R www-data *
sudo chgrp -R development *

./stop.sh

if [ "$arg" == "develop" ]
then
    loglevel="DEBUG"
    config="development.ini"
else
    loglevel="INFO"
    config="production.ini"
fi

echo "Starting workers..."
/data/pyramid/bin/pceleryd $config -Q celery --concurrency=11 --logfile=logs/ptmworker.log --loglevel=$loglevel --pidfile=logs/ptmworker.0.pid > logs/ptmworker.0.start &
/data/pyramid/bin/pceleryd $config -Q notify --concurrency=1 --logfile=logs/ptmworker.log --loglevel=$loglevel --pidfile=logs/ptmworker.1.pid > logs/ptmworker.1.start &

echo "Restarting webservice..."
sudo apache2ctl restart
