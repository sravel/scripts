#! /usr/bin/python
# -*- coding: utf-8 -*-
# Python v3
 
import sys
import os
from PyQt5 import (QtWidgets, QtCore)
 
class Arbodisk(QtWidgets.QWidget):
 
    def __init__(self):
        super().__init__()
        self.resize(400, 300)
 
        # répertoire racine pour le QTreeView
        rephome = QtCore.QDir.homePath()
 
        #-- Modèle
        self.myModel = QtWidgets.QFileSystemModel()
        self.myModel.setReadOnly(False)
        self.myModel.setRootPath(rephome)
 
        #-- treeview
        self.myTreeView = QtWidgets.QTreeView(self)
        self.myTreeView.setModel(self.myModel)
 
        rootModelIndex = self.myModel.index(rephome)
        self.myTreeView.setRootIndex(rootModelIndex)
 
        # positionnement dans le fenêtre
        posit = QtWidgets.QGridLayout()
        posit.addWidget(self.myTreeView, 0, 0)
        self.setLayout(posit)
 
##############################################################################
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    fen = Arbodisk()
    fen.show()
    end = app.exec_()
