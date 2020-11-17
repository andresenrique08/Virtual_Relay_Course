from Curvas import TC_curve as curvaRele


class Function51:
    def __init__(self, T_spam, I_pick, Td, CurveDelay, Index=1, Alpha=0.02, Beta=0.14, L=0, C=0, *name):
        self.__curve = curvaRele.TCCurve(I_pick, Td, Index, Alpha, Beta, L, C, name)
        # curve, is a TC_curve.py object
        self.__timeEl = 0  # time elapsed in seconds
        self.__sampleTime = T_spam  # sample time in seconds
        self.__iCut = I_pick  # Pick up current in Amperes
        self.__delayCurve = CurveDelay  # seconds

    def __addTime(self):
        self.__timeEl = self.__timeEl + self.__sampleTime

    def __restartCounter(self):
        self.__timeEl = 0

    def __checkCurrent(self, I_measured):
        if abs(1000*I_measured) > self.__iCut:
            self.__addTime()
        else:
            self.__restartCounter()

    def getStatus(self, I_measured):
        #####List of messages####
        # return 0, if 'Normal' state
        # return 1, if 'Stand By' state
        # return 2, if 'Lock out' state
        msg = 0
        self.__checkCurrent(I_measured)
        if self.__iCut < 1000 * abs(I_measured):
            tCut = self.__curve.mathCurve(1000 * abs(I_measured)) + self.__delayCurve  # miliseconds
            if self.__timeEl >= tCut:
                msg = 2
            else:
                msg = 1
        return msg
