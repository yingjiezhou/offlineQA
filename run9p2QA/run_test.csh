#!/bin/bash

starver SL19b
cons -r; cons

if [ ! -d test ];then
   mkdir test 
fi

rm -f test/runtestLog
rm -f test/test.root
root4star -b -l <<EOF >& test/runtestLog
.O2
.x doEvent.C(-1, "test.list","test/test.root", 0) 
.q
EOF
