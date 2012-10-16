/data/pyramid/bin/python setup.py build
/data/pyramid/bin/python setup.py develop
sudo chmod -R 775 *
sudo chown -R www-data *
sudo chgrp -R development *

kill `cat "logs/pids"`
/data/pyramid/bin/pceleryd development.ini --logfile=logs/pceleryd --loglevel=DEBUG --pidfile=logs/pids > /dev/null &

sudo apache2ctl restart

