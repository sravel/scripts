#!/bin/bash
# -*- coding: utf-8 -*-
## @package epydocUpdate.sh
# @author Sebastien Ravel


cp ./modules/DOC/index.html ./

epydoc --html ./modules/MODULES_SEB.py -o ./modules/DOC/ -v

mv ./index.html ./modules/DOC/

git commit -m "update DOC for MODULE_SEB" ./modules/DOC/*

scp ./modules/DOC/* admin@194.254.138.95:/volume1/web/DOC/
