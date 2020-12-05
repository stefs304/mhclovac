from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import StackingRegressor, StackingClassifier
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier
from sklearn.svm import SVR, SVC
from .preprocessing import get_features, transform_ic50_measures


class BindingModel:

    model = None
    index_id_list = ['PRAM900101', 'FASG760101', 'ZIMJ680104', 'CHOP780201', 'PRAM820103', 'RACS820112', 'ROBB760107']

    def __init__(self, n_jobs=-1, index_id_list=None, verbose=False):
        if index_id_list:
            self.index_id_list = index_id_list
        self.init_model(n_jobs=n_jobs, verbose=verbose)

    def init_model(self, n_jobs, verbose):
        """
        Create stacked model for sequence features prediction
        """
        level0 = list()
        level0.append(('knn', KNeighborsRegressor()))
        level0.append(('dt', DecisionTreeRegressor()))
        level0.append(('svm', SVR()))
        level1 = LinearRegression()
        stacked_model = StackingRegressor(estimators=level0, final_estimator=level1, cv=5, n_jobs=n_jobs, verbose=verbose)
        self.model = make_pipeline(StandardScaler(), stacked_model)

    def fit(self, peptide_list, ic50_values):
        X = get_features(peptide_list=peptide_list, index_id_list=self.index_id_list)
        y = transform_ic50_measures(ic50_values=ic50_values)
        self.model.fit(X, y)

    def predict(self, peptide_list):
        X = get_features(peptide_list=peptide_list, index_id_list=self.index_id_list)
        return self.model.predict(X)


class LigandModel:

    model = None

    def __init__(self, random_state=None, n_jobs=-1, standardize_data=True):
        pass

    def fit(self, X, y):
        pass

    def predict(self, X):
        pass

    def predict_proba(self, X):
        # if self.scaler:
        #     X = self.scale_transform(X)
        # return [x[1] for x in self.vote_reg.predict_proba(X)]
        pass

    def score(self, X, y):
        pass
