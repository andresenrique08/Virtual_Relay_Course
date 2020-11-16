import numpy as np
import matplotlib.pyplot as plt


class TCCurve:
    def __init__(self, Ip, Td, Indice=1, Alpha=0.02, Beta=0.14, L=0, C=0, *name):
        # indice: indica el método de la curva  | Ip, Pickup current
        # '1', Standard Inverse (IEC)           |
        # '2', Very Inverse (IEC)               |
        # '3', Extremely Inverse (IEC)          |
        # '4', Long Time Standard (IEC)         |
        # '5', Rectifier (UK)                   |
        # '6', Moderately Inverse (IEEE)        |
        # '7', Very Inverse (IEEE)              |
        # '8', Extremely Inverse (IEEE)         |
        # '9', Inverse (US-C08)                 |
        # '10', Short Time Inverse              |
        # '11', Define by user                  |

        if Indice <= 11:
            method_name = 'curve_' + str(Indice)
            method = getattr(self, method_name, lambda: 'Invalid')
            [a, b] = method()
            self.__curveName = a
            self.__alpha = b[0]
            self.__beta = b[1]
            self.__l = b[2]
            self.__c = b[3]

        else:
            self.__curveName = name
            self.__beta = Beta
            self.__alpha = Alpha
            self.__l = L
            self.__c = C
        self.__ip = Ip  # Pickup current
        self.__td = Td  # Time

    # name, beta, alpha, l, c
    def curve_1(self):  # '1', Standard Inverse (IEC)
        return ['Standard Inverse (IEC)', [0.14, 0.02, 0, 0]]

    def curve_2(self):  # '2', Very Inverse (IEC)
        return ['Very Inverse (IEC)', [13.5, 1, 0, 0]]

    def curve_3(self):  # '3', Extremely Inverse (IEC)
        return ['Extremely Inverse (IEC)', [80, 2, 0, 0]]

    def curve_4(self):  # '4', Long Time Standard (IEC)
        return ['Long Time Standard (IEC)', [120, 1, 0, 0]]

    def curve_5(self):  # '5', Rectifier (UK)
        return ['Rectifier (UK)', [45900, 5.6, 0, 0]]

    def curve_6(self):  # '6', Moderately Inverse (IEEE)
        return ['Moderately Inverse (IEEE)', [0.0515, 0.02, 0.114, 0]]

    def curve_7(self):  # '7', Very Inverse (IEEE)
        return ['Very Inverse (IEEE)', [19.61, 2, 0.491, 0]]

    def curve_8(self):  # '8', Extremely Inverse (IEEE)
        return ['Extremely Inverse (IEEE)', [28.2, 2, 0.1217, 0]]

    def curve_9(self):  # '9', Inverse (US-C08)
        return ['Inverse (US-C08)', [5.95, 2, 0.18, 0]]

    def curve_10(self):  # '10', Short Time Inverse
        return ['Short Time Inverse', [0.16758, 0.02, 0.11858, 0]]

    def mathCurve(self, I_measured):
        return (self.__td * (self.__beta / (
                    np.power((I_measured / self.__ip), self.__alpha) - 1) + self.__l) + self.__c) * 1000  # ms

    def graficar(self):
        I = np.linspace(1, 1000, 1000)
        T = self.mathCurve(I)
        print(np.size(T))
        # plt.figure()
        # plt.plot(I, 20*np.log10(T))
        # Add a legend
        plt.loglog(I, T)
        plt.grid(True, which="both", ls="-")
        plt.legend('TC curve')
        plt.xlabel('I')
        plt.ylabel('T')
        plt.title('Curva TC')

        # Show the plot
