import sys

from PyQt4 import uic, QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy as np

from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
import fluorescence as flu


paras = Parameters()
paras.add('a',value=1,min=0,max=2,vary=False)
paras.add('b',value=10,min=8,max=12,vary=True)

pname_to_fit = 'a'
para_to_fit = paras[pname_to_fit]

(ui_err1, QDialog) = uic.loadUiType('err1.ui')
class ErrDialogClass(QDialog,ui_err1):
    def __init__(self, parameters, pname_to_fit, num_interval=20, parent=None):
        QDialog.__init__(self, parent)
        # self.ui = ui_err1()
        self.setupUi(self)
        self.parameters = parameters
        self.para = self.parameters[pname_to_fit]
        self.num_interval = num_interval
        self.init_ui()

        
    def init_ui(self):
        self.bestvalLE.setText(format(self.para.value,'.2e'))
        self.leftLimitLE.setText(format(self.para.value*0.9, '.2e')) # left range: 0.1*value
        self.rightLimitLE.setText(format(self.para.value*1.1, '.2e')) # right range: 0.1*value
        self.numIntervalLE.setText(format(self.num_interval+1,'d')) # number of intervals for the sample


    @pyqtSignature("")
    def fit_range(self):
        best_value, left_limit, right_limit, num_samples = tuple(
            float(self.bestvalLE.text()),
            float(self.leftLimitLE.text()),
            float(self.rightLimitLE.text()),
            float(self.num_samples.text())
        )
        self.fit_range = np.linspace(lift_limit,right_limit,num_samples)
            

    @pyqtSignature("")
    def on_cancelPB_clicked(self):
        print "cancelErrCal"


    def chisq_cal(self):
        chisq = self.rightLimitLE.setText(format(1.87, '.2e'))
        print "calculating chisq"
        return chisq

def main():
    app = QtGui.QApplication(sys.argv)

    def errDialog1()
        err1 = ErrDialogClass(paras,'a')
        err1.nextPB.clicked.connect(errDialog2)
        err1.show()
    def errDialog2()
        err2 = ErrDialogFitChoose()
        err2.nextPB.clidked.connect()
    
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()


