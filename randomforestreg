
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import numpy as np

# Read the Excel file into a DataFrame
df = pd.read_excel("F:/NISER internship/try_2/week_4_new.xlsx")

# Select the features (independent variables) and target (dependent variable)
features = df[['Average_z_score_metabolites']]
target = df['Average']

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(features, target, test_size=0.2, random_state=42)

# Create a random forest regression model
rf = RandomForestRegressor(n_estimators=100, random_state=42)

# Fit the model to the training data
rf.fit(X_train, y_train)

# Predict the target variable for the test set
y_pred = rf.predict(X_test)

# Evaluate the model performance
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

# Print the mean squared error and R-squared value
print("Mean Squared Error:", mse)
print("R-squared:", r2)




# Sort the predicted and actual values based on the index
sorted_index = np.argsort(y_test)
y_test_sorted = np.array(y_test)[sorted_index]
y_pred_sorted = np.array(y_pred)[sorted_index]

# Create a line plot of actual vs predicted values
plt.plot(range(len(y_test)), y_test_sorted, label='Actual')
plt.plot(range(len(y_test)), y_pred_sorted, label='Predicted')
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('week_4')
plt.legend()
plt.show()
