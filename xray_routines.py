import numpy as np

# uipara[0] = float(self.ui.fluxenLE.text()) # energy in KeV
# uipara[1] = float(self.ui.flusliLE.text()) # slit size in mm
# uipara[2] = float(self.ui.fludetLE.text()) # detector length
# uipara[3] = float(self.ui.flubetatopLE.text()) # top beta

def fluCalFun(flupara, refpara, uipara, x):
    surden = flupara[0]  # surface density
    qoff = flupara[1]  # q offset
    yscale = flupara[2]  # y scale
    bgcon = flupara[3]  # background constant
    surcur = flupara[5] * 1e10  # surface curvature, in unit of /AA
    conupbk = flupara[4]  # background linear is borrowed for upper phase concentration.
    conbulk = flupara[6]  # bulk concentration
    k0 = 2 * np.pi * uipara[0] / 12.3984  # wave vector
    slit = uipara[1]  # get slits size
    detlen = uipara[2] * 1e7  # get detector length in unit of /AA
    topd = 1 / (uipara[3] * 2 * k0)  # get the absorption length in top phase: len=1/mu=1/(beta*2*k)
    qz = x + qoff

    refpara['q_off'].value = 0  # reset qoffset in the reflectivity data.
    refpara['rho_b'].value = float(self.ui.flurhobotLE.text())  # set electron density for bottom phase
    refpara['rho_t'].value = float(self.ui.flurhotopLE.text())  # set electron density for top phase

    refModel = self.ref2min(refpara, None, None, None, fit=False, rrf=False)



    alpha = qz / 2 / k0  # get incident angle
    fprint = slit / alpha * 1e7  # get the footprint in unit of /AA
    if surcur == 0:  # no surface curvature
        flu = []
        # p_d=[]
        for i in range(len(alpha)):
            z1 = (fprint[i] - detlen) / 2 * alpha[i]
            z2 = (fprint[i] + detlen) / 2 * alpha[i]
            effd, trans = self.frsnllCal(self.flutopdel, self.flutopbet, self.flubotdel, self.flubotbeta,
                                         self.flubotmu1, k0, alpha[i])
            effv = effd * topd * np.exp(-detlen / 2 / topd) * (detlen * effd * np.exp(z2 / alpha[i] / topd) * (
                        np.exp(-z1 / effd) - np.exp(-z2 / effd)) + topd * (np.exp(detlen / topd) - 1) * (z1 - z2)) / (
                               detlen * effd + topd * (z1 - z2))
            int_sur = surden * topd * (np.exp(detlen / 2 / topd) - np.exp(-detlen / 2 / topd))  # surface intensity
            int_bulk = effv * self.avoganum * conbulk * self.fluelepara[0][
                1] / 1e27  # bluk intensity; the element in the first row is the target element
            int_tot = yscale * trans * (int_sur + int_bulk) + bgcon + bglin * alpha[i]
            flu.append(int_tot)
        # p_d.append(effd)
    else:  # with surface curvature
        flu = []
        self.flu_oil = []
        self.flu_bulk = []
        self.flu_sur = []
        for i in range(len(alpha)):
            bsum = 0
            ssum = 0
            usum = 0
            steps = int((detlen + fprint[i]) / 2 / 1e6)  # use 0.1 mm as the step size
            stepsize = (detlen + fprint[i]) / 2 / steps
            x = np.linspace(-fprint[i] / 2, detlen / 2,
                            steps)  # get the position fo single ray hitting the surface relative to the center of detector area with the step size "steps"
            for j in range(len(x)):
                alphanew = alpha[i] - x[j] / surcur  # the incident angle at position x[j]
                y1 = -detlen / 2 - x[j]
                y2 = detlen / 2 - x[j]
                absorb_y1 = np.exp(-y1 / topd)  # y1 = (-detlen/2-x[j]) distance between x' and left edge of detector
                absorb_y2 = np.exp(-y2 / topd)  # y2 = (detlen/2-x[j]) distance between x' and right edge of detector
                effd, trans = self.frsnllCal(self.flutopdel, self.flutopbet, self.flubotdel, self.flubotbeta,
                                             self.flubotmu1, k0, alphanew)
                ref = refModel(2 * k0 * alphanew)[0]  # calculate the reflectivity at incident angle alpha'.
                absorb_top = np.exp(-x[j] / topd)
                if x[j] > -detlen / 2:
                    bsum = bsum + absorb_top * trans * effd * (1.0 - np.exp(-y2 * alpha[i] / effd))  # equation (5)(1)
                    ssum = ssum + absorb_top * trans
                    # usum = usum + alpha[i] * absorb_top * ((fprint[i]/2-x[j]) + (fprint[i]/2+x[j])*ref[i])
                    usum = usum + absorb_top * alpha[i] * topd * (
                                (1 - 1 / absorb_y1) + ref * (1 - absorb_y2))  # eq (x)(1)
                else:
                    bsum = bsum + absorb_top * trans * effd * (alpha[i] / absorb_y1 / effd) - np.exp(
                        -y2 * alpha[i] / effd))  # #surface has no contribution at this region, equatoin (5)(2)
                    usum = usum + absorb_top * alpha[i] * topd * ref * (
                                np.exp(-y1 / topd) - np.exp(-y2 / topd))  # eq (x)(2)
                    int_bulk = bsum * stepsize * self.avoganum * conbulk * self.fluelepara[0][1] / 1e27
                    int_upbk = usum * stepsize * self.avoganum * conupbk * self.fluelepara[0][
                        1] / 1e27  # if there is metal ions in the upper phase.
                    int_sur = ssum * stepsize * surden
                    int_tot = yscale * (int_bulk + int_sur + int_upbk) + bgcon
                    # int_tot = yscale * (int_bulk + int_sur) + bgcon
                    flu.append(int_tot)
                    self.flu_oil.append(yscale * int_upbk)
                    self.flu_sur.append(yscale * int_sur)
                    self.flu_bulk.append(yscale * int_bulk)
    return flu



def ref2min(params, x, y, yerr, fit=True, rrf=True):
    # residuel for ref fitting
    row = self.ui.refparTW.rowCount()
    d = [params[self.refparaname[i * 4 + 3]].value for i in range(row - 2)]
    rho = [params[self.refparaname[i * 4]].value for i in range(row - 1)]
    mu = [params[self.refparaname[i * 4 + 1]].value for i in range(row - 1)]
    sigma = [params[self.refparaname[i * 4 + 2]].value for i in range(row - 1)]
    rho.append(params[self.refparaname[-2]].value)  # add bottom phase
    mu.append(params[self.refparaname[-1]].value)  # add bottom phase
    syspara = [params[self.refsysparaname[i]].value for i in range(3)]

    if rrf == True:  # whether it is a rrf or ref model
        model = lambda xx: mfit.refCalFun(d, rho, mu, sigma, syspara, xx)
    else:
        model = lambda xx: mfit.refCalFun(d, rho, mu, sigma, syspara, xx, rrf=False)

    if fit == True:  # wether it returns the model or the rsiduals.
        return (model(x) - y) / yerr
    else:
        return model

