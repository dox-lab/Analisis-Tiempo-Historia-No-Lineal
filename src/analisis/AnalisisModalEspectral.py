import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class AnalisisModalEspectral:
    """
    Clase para realizar el análisis dinámico modal espectral en una estructura.
    """
    def __init__(self, zucs_r, matriz_masas, modos, T, num_pisos, num_modos, R, corte_estatico_x, corte_estatico_y, alt_entre_piso):
        """
        Inicializa los parámetros necesarios para el análisis dinámico modal espectral.
        :param zucs_r: Espectro de diseño.
        :param matriz_masas: Matriz de masas del sistema estructural.
        :param modos: Matriz de modos de vibración.
        :param T: Periodos de los modos de vibración.
        :param num_pisos: Número de pisos o niveles de la estructura.
        :param num_modos: Número de modos de vibración a considerar.
        :param R: Factor de reducción sísmica.
        :param corte_estatico_x: Cortante acumulada en la dirección X.
        :param corte_estatico_y: Cortante acumulada en la dirección Y.
        :param alt_entre_piso: Altura entre pisos o niveles de la estructura.
        """
        self.zucs_r = zucs_r
        self.matriz_masas = matriz_masas
        self.modos = modos
        self.T = T
        self.num_pisos = num_pisos
        self.num_modos = num_modos
        self.R = R
        self.corte_estatico_x = corte_estatico_x
        self.corte_estatico_y = corte_estatico_y
        self.alt_entre_piso = alt_entre_piso
        self.total_dof = 3 * num_pisos  # Grados de libertad totales del sistema

    def combinacion_modal(self):
        """
        Realiza la superposición modal espectral y devuelve los resultados en X y Y.
        """
        # Inicialización de las respuestas en X y Y
        disp_abs_x, disp_rmsc_x = np.zeros(self.total_dof), np.zeros(self.total_dof)
        drift_abs_x, drift_rmsc_x = np.zeros(self.total_dof), np.zeros(self.total_dof)
        shear_abs_x, shear_rmsc_x = np.zeros(self.total_dof), np.zeros(self.total_dof)
        disp_abs_y, disp_rmsc_y = np.zeros(self.total_dof), np.zeros(self.total_dof)
        drift_abs_y, drift_rmsc_y = np.zeros(self.total_dof), np.zeros(self.total_dof)
        shear_abs_y, shear_rmsc_y = np.zeros(self.total_dof), np.zeros(self.total_dof)

        # Vectores unitarios en X, Y, y rotación Z
        unit_x, unit_y, unit_rz = np.zeros(self.total_dof), np.zeros(self.total_dof), np.zeros(self.total_dof)
        unit_x[0::3], unit_y[1::3], unit_rz[2::3] = 1, 1, 1

        # Cálculo de las respuestas modales espectrales
        for j in range(len(self.modos)):
            force_x = self.modos[j].T @ self.matriz_masas @ unit_x
            force_y = self.modos[j].T @ self.matriz_masas @ unit_y

            # Cálculo de aceleraciones y desplazamientos espectrales
            spectral_accel = self.zucs_r[j] * 9.80665
            spectral_disp = spectral_accel / (2 * np.pi / self.T[j])**2

            # Respuesta en X
            resp_disp_x = spectral_disp * force_x * self.modos[j]
            resp_accel_x = spectral_accel * force_x * self.matriz_masas @ self.modos[j]
            disp_abs_x += abs(resp_disp_x)
            disp_rmsc_x += resp_disp_x**2
            resp_disp_x[3:] -= resp_disp_x[:-3]
            drift_abs_x += abs(resp_disp_x)
            drift_rmsc_x += resp_disp_x**2
            shear_abs_x += abs(np.cumsum(resp_accel_x[::-1])[::-1])
            shear_rmsc_x += (np.cumsum(resp_accel_x[::-1])[::-1])**2

            # Respuesta en Y
            resp_disp_y = spectral_disp * force_y * self.modos[j]
            resp_accel_y = spectral_accel * force_y * self.matriz_masas @ self.modos[j]
            disp_abs_y += abs(resp_disp_y)
            disp_rmsc_y += resp_disp_y**2
            resp_disp_y[3:] -= resp_disp_y[:-3]
            drift_abs_y += abs(resp_disp_y)
            drift_rmsc_y += resp_disp_y**2
            shear_abs_y += abs(np.cumsum(resp_accel_y[::-1])[::-1])
            shear_rmsc_y += (np.cumsum(resp_accel_y[::-1])[::-1])**2

        # Combinación 25% ABS + 75% RMS para X y Y
        disp_comb_x = 0.25 * disp_abs_x + 0.75 * np.sqrt(disp_rmsc_x)
        drift_comb_x = 0.25 * drift_abs_x + 0.75 * np.sqrt(drift_rmsc_x)
        shear_comb_x = 0.25 * shear_abs_x + 0.75 * np.sqrt(shear_rmsc_x)

        disp_comb_y = 0.25 * disp_abs_y + 0.75 * np.sqrt(disp_rmsc_y)
        drift_comb_y = 0.25 * drift_abs_y + 0.75 * np.sqrt(drift_rmsc_y)
        shear_comb_y = 0.25 * shear_abs_y + 0.75 * np.sqrt(shear_rmsc_y)

        # Crear DataFrame con resultados
        rows = [{
            'Nivel': i + 1,
            'ShearX(kN)': shear_comb_x[0::3][i] / 1000,
            'ShearY(kN)': shear_comb_y[1::3][i] / 1000,
            'DispX(cm)': disp_comb_x[0::3][i] * 100,
            'DispY(cm)': disp_comb_y[1::3][i] * 100
        } for i in range(int(self.total_dof / 3))]

        a_din_sin_escalar = pd.DataFrame(rows)
        return disp_comb_x, drift_comb_x, shear_comb_x, disp_comb_y, drift_comb_y, shear_comb_y, a_din_sin_escalar

    def escalamiento(self, shear_comb_x, shear_comb_y):
        """
        Escala los resultados del análisis dinámico si es necesario y genera mensajes con detalles.
        """
        # Cortantes dinámicos y estáticos
        corte_dinamico_x = shear_comb_x[0::3][0]
        corte_dinamico_y = shear_comb_y[0::3][0]
        corte_estatico_x = 0.80 * self.corte_estatico_x[0]
        corte_estatico_y = 0.80 * self.corte_estatico_y[0]

        # Inicialización de factores de escala y mensajes
        scale_factor_x = 1.0
        scale_factor_y = 1.0
        message_x = f'\nCortante dinámico X: {corte_dinamico_x:.2f} kN, Cortante estático X: {corte_estatico_x:.2f} kN, Porcentaje: {100 * corte_dinamico_x / self.corte_estatico_x[0]:.2f}%. No es necesario escalar en X.'
        message_y = f'Cortante dinámico Y: {corte_dinamico_y:.2f} kN, Cortante estático Y: {corte_estatico_y:.2f} kN, Porcentaje: {100 * corte_dinamico_y / self.corte_estatico_y[0]:.2f}%. No es necesario escalar en Y.'

        # Escalamiento en X si es necesario
        if corte_dinamico_x < corte_estatico_x:
            scale_factor_x = corte_estatico_x / corte_dinamico_x
            message_x = (f'Cortante dinámico X: {corte_dinamico_x:.2f} kN, Cortante estático X: {corte_estatico_x:.2f} kN, '
                        f'Porcentaje: {100 * corte_dinamico_x / self.corte_estatico_x[0]:.2f}%. '
                        f'Factor de escala aplicado en X: {scale_factor_x:.4f}.')

        # Escalamiento en Y si es necesario
        if corte_dinamico_y < corte_estatico_y:
            scale_factor_y = corte_estatico_y / corte_dinamico_y
            message_y = (f'Cortante dinámico Y: {corte_dinamico_y:.2f} kN, Cortante estático Y: {corte_estatico_y:.2f} kN, '
                        f'Porcentaje: {100 * corte_dinamico_y / self.corte_estatico_y[0]:.2f}%. '
                        f'Factor de escala aplicado en Y: {scale_factor_y:.4f}.')

        return scale_factor_x, scale_factor_y, message_x, message_y

    def aplicar_escala(self, scale_factor_x, scale_factor_y, disp_comb_x, drift_comb_x, shear_comb_x, disp_comb_y, drift_comb_y, shear_comb_y):
        """
        Aplica los factores de escala y devuelve los resultados finales.
        """
        rows = [{
            'Nivel': i + 1,
            'F.CortanteX(kN)': scale_factor_x * shear_comb_x[0::3][i] / 1000,
            'F.CortanteY(kN)': scale_factor_y * shear_comb_y[1::3][i] / 1000,
            'DesplazamientoX(cm)': 0.75 * self.R * disp_comb_x[0::3][i] * 100,
            'DesplazamientoY(cm)': 0.75 * self.R * disp_comb_y[1::3][i] * 100,
            'DriftX(‰)': 0.75 * self.R * drift_comb_x[0::3][i] * 1000 / self.alt_entre_piso,
            'DriftY(‰)': 0.75 * self.R * drift_comb_y[1::3][i] * 1000 / self.alt_entre_piso
        } for i in range(self.num_pisos)]

        Adme = pd.DataFrame(rows)
        print('\nANÁLISIS DINÁMICO MODAL ESPECTRAL')
        print(Adme.round(4).to_string(index=False))
        return Adme

    def plot_derivas(self, a_din_sin_escalar):
        """
        Plotea las distorsiones (drifts) en X y Y.
        """
        drift_x = np.array(a_din_sin_escalar['DriftX(‰)'])
        drift_y = np.array(a_din_sin_escalar['DriftY(‰)'])
        limit = 1.1 * max(drift_x.max(), drift_y.max())

        plt.plot(np.insert(drift_x, 0, 0), np.arange(self.num_pisos + 1), 'bo--', label='Drift X', lw=0.8)
        plt.plot(np.insert(drift_y, 0, 0), np.arange(self.num_pisos + 1), 'ro--', label='Drift Y', lw=0.8)
        plt.legend()
        plt.xlabel('Distorsión (‰)')
        plt.ylabel('Nivel')
        plt.axis([-0.05, limit, -0.05, self.num_pisos + 0.05])
        plt.yticks(np.arange(0, self.num_pisos + 0.05, 1))
        plt.show()
