"""
Run and play simulation animation in 2D (need PySide)
"""

import logging
import sys

from PySide import QtGui, QtCore

from kt_simul.gui.game import InteractiveCellWidget

class Animator:
    """
    """

    def __init__(self, meta_instance):
        """
        """
        self.meta = meta_instance

    def play(self,):
        """
        """

        logging.info("Playing animation")

        app = QtGui.QApplication([])
        QtCore.qsrand(QtCore.QTime(0, 0, 0).secsTo(QtCore.QTime.currentTime()))

        widget = InteractiveCellWidget(self.meta)
        widget.show()

        app.exec_()
