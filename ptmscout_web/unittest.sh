if [ $# -eq 0 ]
then
    exitcode=0
    for f in `find tests | grep "test_.*\.py$"`
    do
        proceed=`head -n1 $f | grep '# Skip' | wc -l`
        package=`echo $f | sed "s/\//./g" | sed "s/\(.*\)\.py/\1/"`
        if [ "$proceed" -eq "0" ]
        then
            echo "Running tests: $package"
            PYTHONPATH=/data/ptmscout/ptmscout_web /data/pyramid/bin/python -m unittest $package 2> tmp_test_output.txt
            rc=$?
            output=`cat tmp_test_output.txt`
            rm tmp_test_output.txt
            if [[ $rc != 0 ]]
            then
                echo "Failed test..."
                echo "$output"
                exitcode=$rc
            fi
        else
            echo "Skipping $package"
        fi
    done

    if [ "$exitcode" -eq "0" ]
    then
        echo "Success! All tests passed."
    else
        echo "Failed. Some tests failed."
    fi
    exit $exitcode
else
    PYTHONPATH=/data/ptmscout/ptmscout_web /data/pyramid/bin/python -m unittest $@
fi
