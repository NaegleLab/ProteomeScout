file=$1
shift

/data/pyramid/bin/behave tests/behave/$file --logging-clear-handlers --tags ~@pending $@
