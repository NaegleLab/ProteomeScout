/data/pyramid/bin/python setup.py build > logs/build.log
/data/pyramid/bin/python setup.py develop > logs/deploy.log
sudo chmod -R 775 *
sudo chown -R www-data *
sudo chgrp -R development *

if [ -e "logs/ptmworker.pid" ]
then
    kill `cat "logs/ptmworker.pid"`
fi

/data/pyramid/bin/pceleryd development.ini --logfile=logs/ptmworker.log --loglevel=DEBUG --pidfile=logs/ptmworker.pid > logs/ptmworker.start &
sudo apache2ctl restart

