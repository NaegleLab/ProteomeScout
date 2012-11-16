file=$1
shift

/data/pyramid/bin/behave tests/behave/$file --tags ~@pending $@
