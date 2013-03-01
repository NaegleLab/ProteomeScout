#!/bin/bash

set -e

defmode='install'
arg=${1:-$defmode}

if [ ! -e logs ]
then
    mkdir logs/
fi

echo "Deploying codebase..."
#sudo /data/pyramid/bin/python setup.py build > logs/build.log
sudo /data/pyramid/bin/python setup.py $arg > logs/deploy.log
sudo chmod -R 775 *
sudo chown -R www-data *
sudo chgrp -R development *

echo "Restarting webservice..."
sudo apache2ctl restart
