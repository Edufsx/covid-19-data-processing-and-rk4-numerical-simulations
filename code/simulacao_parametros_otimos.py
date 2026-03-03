import pandas as pd
import matplotlib.pyplot as plt

def S_dot(t, S, I, R, alpha, beta):
    return -alpha*I

def I_dot(t, S, I, R, alpha, beta):
    return (alpha - beta)*I

def R_dot(t, S, I, R, alpha, beta):
    return beta*I
  
def rk4_passo_3por3(t, S_dot, I_dot, R_dot, S, I, R, alpha, beta, h):
    
    # PVI :: S' = S_dot , I' = I_dot e R' = R_dot
    # S_dot, I_dot e R_dot sao as três funcoes que definem o sistema
    # t eh o tempo no instante anterior
    # S, I, R sao as solucoes no instante de tempo anterior
    # h eh o passo de integracao
    
    k1_S = S_dot(t,       S,            I,            R,            alpha, beta)
    k1_I = I_dot(t,       S,            I,            R,            alpha, beta)
    k1_R = R_dot(t,       S,            I,            R,            alpha, beta)
    
    k2_S = S_dot(t + h/2, S + h/2*k1_S, I + h/2*k1_I, R + h/2*k1_R, alpha, beta)
    k2_I = I_dot(t + h/2, S + h/2*k1_S, I + h/2*k1_I, R + h/2*k1_R, alpha, beta)  
    k2_R = R_dot(t + h/2, S + h/2*k1_S, I + h/2*k1_I, R + h/2*k1_R, alpha, beta)     
    
    k3_S = S_dot(t + h/2, S + h/2*k2_S, I + h/2*k2_I, R + h/2*k2_R, alpha, beta)
    k3_I = I_dot(t + h/2, S + h/2*k2_S, I + h/2*k2_I, R + h/2*k2_R, alpha, beta)  
    k3_R = R_dot(t + h/2, S + h/2*k2_S, I + h/2*k2_I, R + h/2*k2_R, alpha, beta)     

    k4_S = S_dot(t + h,   S + h*k3_S,   I + h*k3_I,   R + h*k3_R,   alpha, beta)      
    k4_I = I_dot(t + h,   S + h*k3_S,   I + h*k3_I,   R + h*k3_R,   alpha, beta)   
    k4_R = R_dot(t + h,   S + h*k3_S,   I + h*k3_I,   R + h*k3_R,   alpha, beta)     

    S_novo = S + h*(k1_S + 2*k2_S + 2*k3_S + k4_S)/6    
    I_novo = I + h*(k1_I + 2*k2_I + 2*k3_I + k4_I)/6 
    R_novo = R + h*(k1_R + 2*k2_R + 2*k3_R + k4_R)/6     

    return S_novo, I_novo, R_novo


def rk4_sistemas3por3(S_dot, I_dot, R_dot, condicao_inicial, a, b, alpha, beta, h):
    
    # S_dot, I_dot e R_dot sao as tres funcoes que definem o sistema
    # condicao inicial (x0, y0, z0) para t = a
    # intevalo de integracao [a, b]
    # passo de integracao h
    
    # Passo 0 : Iniciar o vetor de tempo e os vetores para S, I e R
    t_vetor = [a]
    S0, I0, R0 = condicao_inicial
    S_vetor = [S0]
    I_vetor = [I0]
    R_vetor = [R0]    
    
    # Passo 1 :: Calcular N (numero de repeticoes)
    N = int((b-a)/h)

    # Passo 2 :: Executar N vezes o rk4_passo_3por3
    for i in range(N):
        t_novo = t_vetor[i] + h
        S_novo, I_novo, R_novo = rk4_passo_3por3(t_vetor[i], S_dot, I_dot, R_dot, \
                                                 S_vetor[i], I_vetor[i], R_vetor[i],\
                                                 alpha, beta, h)
        t_vetor.append(t_novo)
        S_vetor.append(S_novo)
        I_vetor.append(I_novo)
        R_vetor.append(R_novo)
        
    return t_vetor, S_vetor, I_vetor, R_vetor

def get_parameters(munic, year):
    dado_real = pd.read_csv(f"../data/filtered_data/{munic}_casos_mm7d_{year}")
    b = dc_munic_e_ano[munic][year]["b"]
    condicao_inicial = dc_munic_e_ano[munic][year]["condicoes_iniciais"]
    alpha = dc_munic_e_ano[munic][year]["alpha"]
    beta = dc_munic_e_ano[munic][year]["beta"]

    return dado_real, b, condicao_inicial, alpha, beta

def plot_casos(munic, year):
    dado_real, b, condicao_inicial, alpha, beta  = get_parameters(munic, year)
    t_vector, S_vector, I_vector, R_vector = rk4_sistemas3por3(S_dot, I_dot, R_dot, \
                                               condicao_inicial, a, b, \
                                               alpha, beta, h)
    plt.figure()
    plt.plot(t_vector, I_vector)
    plt.plot(dado_real)
    plt.xlabel(f'Dias - {year}')
    plt.ylabel(f'Número de casos Confirmados em {munic}')
    plt.title(f"Casos de Covid em {munic} (Oficial e Simulado) ao longo de {year}")
    plt.legend(["Simulado", "Oficial"])
    plt.show()

dc_munic_e_ano = {
    "Bauru" : {
        2020 : {"alpha" : 0.576, 
                "beta" : 0.541, 
                "b" : 272,
                "condicoes_iniciais" :[364223, 2, 0]
                },
        2021 : {"alpha" :  0.541, 
                "beta" : 0.516, 
                "b" : 365,
                "condicoes_iniciais" :[345928, 8, 18289]
                },
        2022 : {"alpha" : 0.360, 
                "beta" : 0.330,
                "b": 365,
                "condicoes_iniciais" : [310282, 1, 53942]
                }
},   
    "São Paulo" : {
        2020 : {"alpha" : 0.584, 
                "beta" : 0.541, 
                "b" : 310,
                "condicoes_iniciais" : [11869659, 1, 0]
                },
        2021 : {"alpha" :  0.541, 
                "beta" : 0.522, 
                "b" : 365,
                "condicoes_iniciais" : [11466825, 1117, 401718]
                },
        2022 : {"alpha" : 0.541, 
                "beta" : 0.518,
                "b": 365,
                "condicoes_iniciais" : [10891674, 68, 977918]
                }}
}

# Integracao Numerica
h = 10**(-3)
a = 0

#Plotagem
cidades = ["Bauru", "São Paulo"]
anos = [2020, 2021, 2022]

for i in cidades:
    for k in anos:
        plot_casos(i, k)
