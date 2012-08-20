rm -rf build/*
rm -rf dist/*
rm -rf ptmscout.egg-info/*
/data/pyramid/bin/python setup.py build
/data/pyramid/bin/python setup.py install
sudo chmod -R 775 *
sudo chown -R www-data *
sudo chgrp -R development *
sudo apache2ctl restart
