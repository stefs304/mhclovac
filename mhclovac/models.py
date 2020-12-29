from sklearn.linear_model import LinearRegression, BayesianRidge
from sklearn.ensemble import (
    RandomForestRegressor,
    GradientBoostingRegressor,
    VotingRegressor
)


class BindingModel:

    model = None

    def __init__(self, n_jobs=-1, verbose=False, random_state=None):
        estimators = [
            ('rf', RandomForestRegressor(max_depth=3, random_state=random_state)),
            ('lr', LinearRegression()),
            ('br', BayesianRidge()),
            ('gb', GradientBoostingRegressor(max_depth=4, random_state=random_state))
        ]
        self.estimators = estimators
        self.model = VotingRegressor(estimators, n_jobs=n_jobs, verbose=verbose)

    def fit(self, X, y):
        self.model.fit(X, y)

    def predict(self, X):
        return self.model.predict(X)

