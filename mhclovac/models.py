from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import StackingRegressor, StackingClassifier
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier
from sklearn.svm import SVR, SVC


class BindingModel:

    model = None

    def __init__(self, n_jobs=-1, verbose=False):
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
        stacked_model = StackingRegressor(
            estimators=level0,
            final_estimator=level1,
            cv=5,
            n_jobs=n_jobs,
            verbose=verbose
        )
        self.model = make_pipeline(StandardScaler(), stacked_model)

    def fit(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)


class LigandModel:

    model = None

    def __init__(self, n_jobs=-1, verbose=False):
        self.init_model(n_jobs=n_jobs, verbose=verbose)

    def init_model(self, n_jobs, verbose):
        level0 = list()
        level0.append(('lr', LogisticRegression()))
        level0.append(('knn', KNeighborsClassifier()))
        level0.append(('cart', DecisionTreeClassifier()))
        level0.append(('svm', SVC()))
        level0.append(('bayes', GaussianNB()))
        level1 = LogisticRegression()
        stacked_model = StackingClassifier(
            estimators=level0,
            final_estimator=level1,
            cv=5,
            n_jobs=n_jobs,
            verbose=verbose
        )
        self.model = make_pipeline(StandardScaler(), stacked_model)

    def fit(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)

    def predict_proba(self, X):
        # if self.scaler:
        #     X = self.scale_transform(X)
        # return [x[1] for x in self.vote_reg.predict_proba(X)]
        pass
