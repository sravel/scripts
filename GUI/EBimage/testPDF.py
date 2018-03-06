#!/usr/bin/python3.5
# -*- coding: utf-8 -*-
### @package GUI_EBimage.py
## @author Sebastien Ravel


#import sys

#from PyQt5.QtCore import *
#from PyQt5.QtGui import *
#import popplerqt5
#from pictureflow import *


#w = PictureFlow()
#d = popplerqt5.Poppler.Document.load('file.pdf')
#d.setRenderHint(popplerqt5.Poppler.Document.Antialiasing and popplerqt5.Poppler.Document.TextAntialiasing)

#page = 0
#pages = d.numPages() - 1
#while page < pages:
	#page += 1
	#print(page)
	#w.addSlide(d.page(page).renderToImage())
#w.show()

#sys.exit(app.exec_())
# -*- coding: utf-8 -*-

import sys
from PyQt5 import QtGui, QtCore
from truc import Ui_fenetre_principale,PreviewWindow

class MonAppli(QtGui.QMainWindow, Ui_fenetre_principale):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_fenetre_principale.__init__(self)

        # Configure l'interface utilisateur.
        self.setupUi(self)
        self.connect(self.generer, QtCore.SIGNAL("clicked()"), self.showPrev)

    def showPrev(self):
        Dialog = QtGui.QDialog(self)
        ui = PreviewWindow(Dialog)
        Dialog.show()

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MonAppli()
    window.show()
    sys.exit(app.exec_())
