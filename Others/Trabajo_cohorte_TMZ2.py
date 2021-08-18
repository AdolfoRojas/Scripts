#!/usr/bin/env python3
from operator import le
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import sklearn.preprocessing
import sklearn.model_selection
import sklearn.ensemble
import sklearn.linear_model
import sklearn.neighbors
import sklearn.pipeline # Módulo para el establecimento fácil de flujos de trabajo
import sklearn.preprocessing # Preparación de los conjuntos de datos para entrenamiento
import sklearn.model_selection # Funciones para validar métodos
import sklearn.dummy
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut
import plotly.express as px # Biblioteca de visualización

data = pd.read_csv("Otros_proyectos/Cohorte_TMZ/ML_TZM.tab", sep = '\t').drop(columns=['ID','X6mwt'])
#data = pd.read_csv("https://raw.githubusercontent.com/fvillena/biocompu/master/data/insurance.csv")
data.describe()
correlation = data.corr(method="spearman") # Por defecto este método calcula la correlación de Pearson.
correlation
sns.set(font_scale=0.3)
corr_plot = sns.heatmap(correlation, cmap="RdYlGn", center = 0,xticklabels = 1, yticklabels =1) # RdYlGn Utilizamos una paleta de colores divergente para visualizar diferencias en los valores positivos y negativos
corr_plot.get_figure().savefig("Otros_proyectos/Cohorte_TMZ/HeatMap_spearman_TMZ.png", dpi=300)

features = data.iloc[:,:-1]
feature_names = features.columns
label = data["TMZ_6mwt_treat_eff"]
#label = data["charges"]

scaler = sklearn.preprocessing.MinMaxScaler() # Instanciamos la clase del scaler que realizará una normalización de los datos
features[:] = scaler.fit_transform(features) # Ajustamos el scaler y al mismo tiempo transformamos nuestras características

########################################################################################################################################

features_train, features_test, label_train, label_test = sklearn.model_selection.train_test_split( 
    features,
    label,
    test_size=0.30, 
    random_state = 11, 
)

def regression_report(y_true, y_pred):
    """
    Esta función recibe un arreglo de valores reales y predichos para 
    retornar un diccionario con una serie de métricas de regresión
    """
    return {
        'mae': sklearn.metrics.mean_absolute_error(y_true, y_pred),
        'rmse': sklearn.metrics.mean_squared_error(y_true, y_pred) ** 0.5,
        'r2': sklearn.metrics.r2_score(y_true, y_pred)
    }
lr = sklearn.linear_model.LinearRegression()
lr.fit(features_train, label_train)

lr_predictions = lr.predict(features_test)


lr_vil = pd.DataFrame(lr.coef_,columns=["value"]) # Guardamos en un dataframe los coeficientes
lr_vil.index = feature_names
lr_vil

lr_regression_report = regression_report(label_test, lr_predictions)
lr_regression_report

knn = sklearn.neighbors.KNeighborsRegressor() # Instanciamos una support vector machine con un kernel lineal
knn.fit(features_train, label_train)

knn_predictions = knn.predict(features_test)

knn_regression_report = regression_report(label_test, knn_predictions)
knn_regression_report

rf = sklearn.ensemble.RandomForestRegressor()
rf.fit(features_train, label_train)
rf_predictions = rf.predict(features_test)
rf_vil = pd.DataFrame(list(zip(feature_names,rf.feature_importances_)),
             columns=["feature","importance"]
            ).set_index("feature")
rf_vil.sort_values("importance",ascending=False)

rf_regression_report = regression_report(label_test, rf_predictions)
rf_regression_report

dummy = sklearn.dummy.DummyRegressor()
dummy.fit(features_train, label_train)

dummy_predictions = dummy.predict(features_test)

dummy_regression_report = regression_report(label_test, dummy_predictions)
dummy_regression_report

performances = pd.DataFrame( # Consolidamos todas las métricas en un DatFrame
    data = [
        lr_regression_report,
        knn_regression_report,
        rf_regression_report,
        dummy_regression_report
    ],
    index = [
        "Linear Regression",
        "k-Nearest Neighbors",
        "Random Forest",
        "Regresor Tonto"
    ]
).sort_values( # Ordenamos los valores
    by="rmse"
)
performances

########################################################################################################################


def make_pipe(regressor):
    """Esta función recibe un estimador para regresión y retorna un Pipeline que antepone un escalador"""
    pipe = sklearn.pipeline.Pipeline( # Establecemos el flujo de trabajo
        [
            ('scaler', sklearn.preprocessing.MinMaxScaler()), # Primero escalamos las características
            ('regressor', regressor) # Después pasamos los datos por el regresor
        ]
    )
    return pipe

regressors = [
    sklearn.linear_model.LinearRegression(), # Regresión Lineal
    sklearn.neighbors.KNeighborsRegressor(), # kNN
    sklearn.ensemble.RandomForestRegressor(random_state = 11) # Random Forest
]


param_grids = [
    { # Regresión Lineal
        'regressor__fit_intercept': [True,False] # En la regresión lineal ajustaremos el intercepto o no.
    },
    { # kNN
        'regressor__n_neighbors': range(2,10,1), # Número de vecinos que se utilizarán para predecir
        'regressor__metric': ['euclidean', 'manhattan'] # Métrica de distancia que se utilizará para encontrar los vecinos más cercanos.
    },
    { # Randos Forest
        'regressor__n_estimators': range(100,1001,300), # Cantidad de árboles que se entrenarán
        'regressor__max_features': ['sqrt', 'log2', None] # Cantidad de características que se utilizarán para entrenar cada arbol.
    }
]


results = []
for regressor, param_grid in zip(regressors, param_grids):
    pipe = make_pipe(regressor)
    gs = sklearn.model_selection.GridSearchCV( # Grid search
        pipe, # Este es el estimador para el cual estamos buscando los mejores hiperparámetros
        param_grid, # Este es el espacio de hiperparámetros donde buscaremos.
        scoring = "r2", # Esta es la métrica de desempeño que optimizaremos
        n_jobs = -1, # Esta es la cantidad de trabajos que se calcularán en paralelo.
        cv = 3, # Esta es la cantidad de particiones que tendrá la validación cruzada de cada combinación de hiperparámetros
        verbose = 3 # Qué tanta información se mostrará en la búsqueda.
    )
    gs.fit(features,label)
    results.append(gs)

for result in results:
    print(result.best_estimator_["regressor"].__class__.__name__, result.best_params_,result.best_score_)



knn_results = pd.DataFrame(results[1].cv_results_).round(4).sort_values("mean_test_score")
px.parallel_categories(
    knn_results.filter(regex=r"param_regressor.*|mean_test_score"),color = "mean_test_score").write_image("Otros_proyectos/Cohorte_TMZ/kNN_optimizacion.png")

rf_results = pd.DataFrame(results[2].cv_results_).round(4).sort_values("mean_test_score")
px.parallel_categories(rf_results.filter(regex=r"param_regressor.*|mean_test_score"),color = "mean_test_score").write_image("Otros_proyectos/Cohorte_TMZ/RF_optimizacion.png")

cv_results = []
for result in results:
    best_estimator = result.best_estimator_
    estimator_name = result.best_estimator_["regressor"].__class__.__name__
    scores = sklearn.model_selection.cross_val_score( # Esta función nos retorna una lista de rendimientos para cada partición.
        best_estimator,
        features,
        label, 
        scoring="r2",
        #scoring="neg_mean_absolute_error",  
        n_jobs=-1, 
        cv = 3,
        #cv = 27,  
        verbose = 3)
    cv_results.extend(zip([estimator_name]*len(scores),range(len(scores)), scores))
cv_results = pd.DataFrame(cv_results, columns = ["estimator_name","fold","score"])


cv_results.sample(10)
cv_results.groupby("estimator_name").agg({"score":["mean","std","min","max"]})


px.box(cv_results,x = "estimator_name",y = "score").write_image("Otros_proyectos/Cohorte_TMZ/boxplot_cv_10_comparativo.png")