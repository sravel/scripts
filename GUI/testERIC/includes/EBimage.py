# -*- coding: utf-8 -*-

"""
Module implementing GUI_EBimage.
"""

from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QMainWindow

from .Ui_EBimage import Ui_GUI_EBimage


class GUI_EBimage(QMainWindow, Ui_GUI_EBimage):
    """
    Class documentation goes here.
    """
    def __init__(self, parent=None):
        """
        Constructor
        
        @param parent reference to the parent widget
        @type QWidget
        """
        super(GUI_EBimage, self).__init__(parent)
        self.setupUi(self)
