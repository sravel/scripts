#!/bin/bash
# -*- coding: utf-8 -*-
## @package copyLocalToCluster.sh
# @author Sebastien Ravel



cp ./local/* ./cluster/


for f in ./cluster/*.py
do

	sed -i -e "s|#!/usr/bin/python3.*|#!/usr/local/bioinfo/python/3.4.3_build2/bin/python|" $f
	sed -i -e "s|#!/usr/bin/python2.*|#!/usr/local/bioinfo/python/2.7.9_build2/bin/python|" $f
done


for f in ./cluster/*.R
do

	sed -i -e "s|#!/usr/bin/Rscript --vanilla|#!/usr/local/bioinfo/R/3.2.2/bin/Rscript --vanilla|" $f
done
