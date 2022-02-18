# COVID19-SEIRD-model-for-Africa-2021
This model analyses data for dynamic modelling of the true burden of COVID-19 and deaths in the African region. 
This being country specific modeling, the burden in each of the countries is computed using the Partially Observed Markov Processes (POMP). 
The software used is pomp King et al., 2016 and the data is also attached herein.

#Get started
1.	Create a working directory, which may be local or online
2.	Remove any object that could be in the memory of the R server
3.	Create a data vector containing each of the 47 countries of interest in Africa. This will be used for looping (from country to country)
4.	Load the packages needed, and set the random seed, to allow reproducibility
5.	Select the countries and load the two datasets: CovCasesWeekly.txt that contains the COVID-19 reported cases and deaths; DynCovInitGA.txt that contains the initialization values for all the parameters used in the model.
6.	Filter the data by each of the countries and fit the model process

#7.	In the modelling phase:

a.	Fit cases and deaths data to obtained the measurement model parameters for the distribution assumed 

b.	Initializing the parameters

c.	Derive immunization rates proportions of the population at any time t, using a logistic growth curve

d.	Obtain the seasonality in the data by fitting a Fourier transform to the observed data. The parameters form this model are used to introduce the waves in the force of 
infection rate at the same points as in the observed data

e.	Fit the measurement model C snippets by defining random simulator of measurement model with the parameters obtained from the data

f.	Initialize and define the process as per the state-space model developed

g.	Simulate and optimize the model based on the initialization parameters provided

h.	Aggregate the total number of persons under either sub-model, that is, immunized and not immunized arms and get the 95% confidence intervals

i.	Codes for plotting the same are also provided

8.	The model loop continues and stops when the list of countries has been exhausted
9.	Projection for 2022 starting from the current immunized populations per country at the current rate of immunization and rates of immunization and variants follows the same modelling cycle but now starts at the populations immunized for each country
10.	case scenarios
11.	Re – infection: A new variant that is more re-infectious than the delta by -80%, 120% and 200% assuming dominance of the said new variant
12.	Variant severity:  A new variant that is more severe than the delta by 80%, 120%, and 200% assuming dominance of the said new variant
13.	Transmissibility:  A new variant that is more transmissible than the delta by 80%, 120%, and 200% assuming dominance of the said new variant

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
