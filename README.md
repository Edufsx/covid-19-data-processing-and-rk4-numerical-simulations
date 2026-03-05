# (EN 🇺🇸) COVID-19 Data Analysis and Numerical Simulations of the SIR Model

## Introduction
This repository contains **Python** codes for processing official COVID-19 data and simulating the pandemic dynamics using the **Susceptible-Infected-Recovered (SIR) epidemiological model**, solved numerically using the **fourth-order Runge–Kutta method (RK4)**.

COVID-19 data from the municipalities of **São Paulo** and **Bauru** during the years **2020, 2021, and 2022** were selected. The official data were obtained from the São Paulo State Data Analysis System (SEADE), and the CSV file containing the raw data required for this code can be found at: [COVID-19 Data for São Paulo from SEADE](https://repositorio.seade.gov.br/dataset/covid-19/resource/d2bad7a1-6c38-4dda-b409-656bff3fa56a).

In the link below, additional formats for the SEADE COVID-19 database can also be found: https://repositorio.seade.gov.br/dataset/covid-19.



## Project Structure
- **data/**: contains the raw data and the filtered data used in the simulations;  
- **data_filter.py**: script responsible for processing and filtering the data;  
- **sir_model_simulation.py**: implementation of the SIR model and the numerical simulations.



## Data Processing

The filtering script performs the following steps:

- Selection of data by **municipality** and **year**;  
- Adjustment of the number of cumulative cases considering the final value from the previous year, thus converting the historical cumulative series into annual cumulative values;  
- Removal of initial values equal to zero;  
- Calculation of the **7-day moving average** of cumulative cases to smooth daily fluctuations;  
- Export of the filtered data to `.csv` files.

These processed data are later used in the epidemiological model simulations.



## Epidemiological Model

The **SIR model** describes the dynamics of an epidemic by dividing the population into three compartments:

- **S (Susceptible)** – individuals who can contract the disease;
- **I (Infected)** – individuals currently infected;  
- **R (Recovered)** – individuals who have recovered from the disease. 

Additionally, the parameters considered in the model are:

- Infection rate ($\alpha$) - determines the speed at which a disease spreads;  
- Recovery rate ($\beta$) - determines the speed at which infected individuals recover from the disease or die before recovering;  
- Rate at which the susceptible population is transferred to the infected population $\left(\dfrac{\alpha IS}{N}\right)$;  
- Rate at which the infected population is transferred to the removed population ($\beta{I}$).

The temporal evolution of these compartments is described by a system of **ordinary differential equations**, solved numerically using the **fourth-order Runge–Kutta method (RK4)**.

The system described is the following:

$$
\begin{cases}
    \dfrac{dS}{dt}=\,\,\,\,\,\dfrac{-\alpha IS}{N} \\\\
    \dfrac{dI}{dt}=\,\,\,\,\,\dfrac{\alpha IS}{N} -\beta{I} \\\\
    \dfrac{dR}{dt}=\,\,\,\,\,\beta{I}
\end{cases}
$$

The initial conditions of the model at time $t_0$ are defined as: $S(t)=S_0$, $I(t)=I_0$ and $R(t)=R_0$; and in the simulations it was considered reasonable to assume $S(t)\approx{N}$ because, according to the authors (FRANCO; DUTRA, 2021), for a constant population $S(t)\leq{S}_0$, the number of infected individuals $I(t)$ and recovered individuals $R(t)$ is small compared to the total population $N$. In this way, the model described above can be rewritten as

$$
\begin{cases}
\dfrac{dS}{dt} = -\alpha I \\[6pt]

\dfrac{dI}{dt} = (\alpha - \beta) I \\[6pt]

\dfrac{dR}{dt} = \beta I
\end{cases}
$$

**Reference**:

FRANCO, C. M. R.; DUTRA, R. F. **SIR model for propagation of COVID-19 in the Paraíba's State (Brazil)**. *Intermaths*, v. 2, n. 2, p. 39-48, 2021.

## Code Features

- Implementation of the **SIR epidemiological model**;
- Numerical integration using the **fourth-order Runge–Kutta method**;
- Search for model parameters (**infection and recovery rates**);
- Estimation of the initial number of infected individuals;
- Graphical comparison between **simulated data and observed data**.



## Libraries Used

- `pandas`
- `numpy`
- `matplotlib`
- `tqdm`


## Notes

Some parameter search routines perform a large number of model simulations and may require **high computational time** to execute.

## Author

**Eduardo Ferreira da Silva**

---

# (PT 🇧🇷) Análise de Dados da COVID-19 e Simulações Numéricas do Modelo SIR 

## Introdução
Este repositório contém códigos em **Python** para processamento de dados oficiais da COVID-19 e simulação da dinâmica pandêmica utilizando o **modelo epidemiológico Suscetível-Infectado-Recuperado (SIR)** resolvido numericamente pelo **método de Runge–Kutta de quarta ordem (RK4)**.

Escolheu-se os dados de COVID-19 dos municípios de **São Paulo** e **Bauru** durante os anos de **2020, 2021 e 2022**. Os dados oficiais foram consultados no Sistema Estadual de Análise de Dados (SEADE) e o arquivo CSV que contém os dados brutos necessários para esse código pode ser encontrado em: [Dados COVID19 para São Paulo do SEADE](https://repositorio.seade.gov.br/dataset/covid-19/resource/d2bad7a1-6c38-4dda-b409-656bff3fa56a).

Já no link abaixo, encontra-se mais formatos para a base de dados de COVID-19 do SEADE: https://repositorio.seade.gov.br/dataset/covid-19.


## Estrutura do Projeto
- **data/**: contém os dados brutos e os dados filtrados utilizados nas simulações; 
- **data_filter.py**: script responsável pelo tratamento e filtragem dos dados;  
- **sir_model_simulation.py**: implementação do modelo SIR e das simulações numéricas.


## Processamento de Dados

O script de filtragem realiza as seguintes etapas:

- Seleção dos dados por **município** e **ano**;  
- Ajuste do número de casos acumulados considerando o valor final do ano anterior, assim passando de acumulado da série histórica para acumulado anual;
- Remoção de valores iniciais iguais a zero  
- Cálculo da **média móvel de 7 dias** dos casos acumulados para suavizar flutuações diárias  
- Exportação dos dados filtrados para arquivos `.csv`

Esses dados tratados são utilizados posteriormente nas simulações do modelo epidemiológico.


## Modelo Epidemiológico

O **modelo SIR** descreve a dinâmica de uma epidemia dividindo a população em três compartimentos:

- **S (Suscetíveis)** – indivíduos que podem contrair a doença;  
- **I (Infectados)** – indivíduos atualmente infectados;  
- **R (Recuperados)** – indivíduos que se recuperaram da doença.  

Além disso, os parâmetros considerados no modelo são:

- Taxa de infecção ($\alpha$) - determina a velocidade de propagação de uma doença; 
- Taxa de recuperação ($\beta$) - determina a velocidade com que indivíduos infectados se recuperam da doença ou morrem antes de se recuperar; 
- Taxa que a população suscetível está sendo transferida para a população infectada $\left(\dfrac{\alpha IS}{N}\right)$;
- Taxa que a população infectada está sendo transferida para a população removida ($\beta{I}$).

A evolução temporal desses compartimentos é descrita por um sistema de **equações diferenciais ordinárias**, resolvido numericamente pelo **método de Runge–Kutta de quarta ordem (RK4)**.

O sistema descrito é o seguinte:
$$
\begin{cases}
    \dfrac{dS}{dt}=\,\,\,\,\,\dfrac{-\alpha IS}{N} \\\\
    \dfrac{dI}{dt}=\,\,\,\,\,\dfrac{\alpha IS}{N} -\beta{I} \\\\
    \dfrac{dR}{dt}=\,\,\,\,\,\beta{I}
\end{cases}
$$

As condições iniciais do modelo no momento $t_0$ são definidas como: $S(t)=S_0$, $I(t)=I_0$ e $R(t)=R_0$; e nas simulações foi considerado razoável assumir $S(t)\approx{N}$, pois, segundo os autores (FRANCO; DUTRA, 2021) para uma população constante $S(t)\leq{S}_0$, o número de infectados $I(t)$ e o número de recuperados $R(t)$ é pequeno comparado com o número do total da população $N$. Desta maneira, o modelo descrito acima pode ser reescrito como

$$
\begin{cases}
\dfrac{dS}{dt} = -\alpha I \\[6pt]

\dfrac{dI}{dt} = (\alpha - \beta) I \\[6pt]

\dfrac{dR}{dt} = \beta I
\end{cases}
$$

**Referência**:

FRANCO, C. M. R.; DUTRA, R. F. **SIR model for propagation of COVID-19 in the Paraíba's State (Brazil)**. *Intermaths*, v. 2, n. 2, p. 39-48, 2021.

## Funcionalidades do Código

- Implementação do modelo epidemiológico **SIR**;
- Integração numérica utilizando **Runge–Kutta de quarta ordem**;
- Busca de parâmetros do modelo (**taxa de infecção e recuperação**);
- Estimativa do número inicial de infectados;
- Comparação gráfica entre **dados simulados e dados observados**.


## Bibliotecas Utilizadas

- `pandas`
- `numpy`
- `matplotlib`
- `tqdm`


## Observações

Algumas rotinas de busca de parâmetros realizam um grande número de simulações do modelo e podem demandar **tempo computacional elevado** para serem executadas.

## Autor

**Eduardo Ferreira da Silva**