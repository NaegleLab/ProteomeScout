if [ $# -eq 0 ]
then
    for f in `find tests | grep "test_.*\.py$"`
    do
        proceed=`head -n1 $f | grep '# Skip' | wc -l`
        if [ "$proceed" -eq "0" ]
        then
            package=`echo $f | sed "s/\//./g" | sed "s/\(.*\)\.py/\1/"`
            echo "Running tests: $package"
            PYTHONPATH=/data/ptmscout/ptmscout_web /data/pyramid/bin/python -m unittest $package 2> tmp_test_output.txt
            rc=$?
            output=`cat tmp_test_output.txt`
            rm tmp_test_output.txt
            if [[ $rc != 0 ]]
            then
                echo "Failed test..."
                echo "$output"
                exit $rc
            fi
        fi
    done

    echo "Success!"
else
    PYTHONPATH=/data/ptmscout/ptmscout_web /data/pyramid/bin/python -m unittest $@
fi
