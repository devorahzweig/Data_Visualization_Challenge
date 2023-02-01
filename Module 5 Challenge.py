#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.stats import linregress


# In[8]:


#Preparing the Mouse Data
mouse_data = pd.read_csv("Mouse_metadata.csv")
mouse_data.head()


# In[20]:


unique_mice = mouse_data["Mouse ID"].unique()
unique_mice


# In[11]:


merged = pd.merge(mouse_data, study_results, how='outer', on="Mouse ID")
merged.head()


# In[12]:


unique_mice = merged["Mouse ID"].unique
unique_mice


# In[15]:


#Remove Duplicates 
duplicate_ID = merged.loc[merged.duplicated(subset=['Mouse ID', 'Timepoint']),'Mouse ID'].unique()
df = merged[merged['Mouse ID'].isin(duplicate_ID)==False]
df.head()


# In[18]:


new_unique = len(df["Mouse ID"].unique())
new_unique


# In[22]:


#Generate Summary Statistics
median = df['Tumor Volume (mm3)'].groupby(df['Drug Regimen']).median()
mean = df['Tumor Volume (mm3)'].groupby(df['Drug Regimen']).mean()
var = df['Tumor Volume (mm3)'].groupby(df['Drug Regimen']).var()
std = df['Tumor Volume (mm3)'].groupby(df['Drug Regimen']).std()
SEM = df['Tumor Volume (mm3)'].groupby(df['Drug Regimen']).sem()
summary = pd.DataFrame({"Median Value":median, 
                        "Mean Value":mean, 
                        "Variance":var, 
                        "Std. Dev":std, 
                        "Std. Err":SEM})
summary.head()


# In[28]:


#Bar Chart #1 Using Pandas
mice_count = df["Drug Regimen"].value_counts()
pandas_bar = mice_count.plot(kind = "bar", color = 'orange')
plt.xlabel("Drug Regimen")
plt.ylabel("Mine Count")
plt.title("Mice per Treatment")


# In[30]:


#Bar Chart #2 Using MatPlotLib
x_values = mice_count.index.values
y_values = mice_count.values
plt.bar(x_values, y_values, color='orange', alpha=0.8, align='center')
plt.title("Mice per Treatment")
plt.xlabel("Drug Regimen")
plt.ylabel("Mice Count")
plt.xticks(rotation="vertical")
plt.show()


# In[31]:


#Pie Chart #1 Using Pandas
gender_data = df["Sex"].value_counts()
gender_data.plot.pie(autopct= "%1.1f%%")
plt.title("Female vs. Male Mice")
plt.show()


# In[32]:


#Pie Chart #2 Using MatPlotLib
labels = ['Female', 'Male']
sizes = [49, 51]
plot = gender_data.plot.pie(y='Total Count', autopct="%1.1f%%")
plt.title('Female vs. Male Mice')
plt.ylabel('Sex')
plt.show()


# In[34]:


Capomulin = df.loc[df["Drug Regimen"] == "Capomulin",:]
Ramicane = df.loc[df["Drug Regimen"] == "Ramicane", :]
Infubinol = df.loc[df["Drug Regimen"] == "Infubinol", :]
Ceftamin = df.loc[df["Drug Regimen"] == "Ceftamin", :]


# In[35]:


#Capomulin Data
Capomulin_last = Capomulin.groupby('Mouse ID').max()['Timepoint']
Capomulin_volume = pd.DataFrame(Capomulin_last)
Capomulin_merge = pd.merge(Capomulin_volume, df, on=("Mouse ID","Timepoint"),how="left")
Capomulin_merge.head()


# In[38]:


Capomulin_tumors = Capomulin_merge["Tumor Volume (mm3)"]
quartiles =Capomulin_tumors.quantile([.25,.5,.75])
lowq = quartiles[0.25]
upq = quartiles[0.75]
iqr = upq-lowq
print(f"The lower quartile of Capomulin tumors is: {lowq}")
print(f"The upper quartile of Capomulin tumors is: {upq}")
print(f"The interquartile range of Capomulin tumors is: {iqr}")
lower_bound = lowq - (1.5*iqr)
upper_bound = upq + (1.5*iqr)
print(f"Numbers less than {lower_bound} might be outliers.")
print(f"Numbers more than {upper_bound} might be outliers.")


# In[39]:


#Ramicane Data
Ramicane_last = Ramicane.groupby('Mouse ID').max()['Timepoint']
Ramicane_vol = pd.DataFrame(Ramicane_last)
Ramicane_merge = pd.merge(Ramicane, df, on=("Mouse ID","Timepoint"),how="left")
Ramicane_merge.head()


# In[41]:


Ramicane_tumors = Ramicane_merge["Tumor Volume (mm3)_x"]
quartiles =Ramicane_tumors.quantile([.25,.5,.75])
lowq = quartiles[0.25]
upq = quartiles[0.75]
iqr = upq-lowq
print(f"The lower quartile of Ramicane tumors is: {lowq}")
print(f"The upper quartile of Ramicane tumors is: {upq}")
print(f"The interquartile range of Ramicane tumors is: {iqr}")
lower_bound = lowq - (1.5*iqr)
upper_bound = upq + (1.5*iqr)
print(f"Numbers less than {lower_bound} might be outliers.")
print(f"Numbers more than {upper_bound} might be outliers.")


# In[42]:


#Infubinol Data
Infubinol_last = Infubinol.groupby('Mouse ID').max()['Timepoint']
Infubinol_vol = pd.DataFrame(Infubinol_last)
Infubinol_merge = pd.merge(Infubinol, df, on=("Mouse ID","Timepoint"),how="left")
Infubinol_merge.head()


# In[43]:


Infubinol_tumors = Infubinol_merge["Tumor Volume (mm3)_x"]
quartiles =Infubinol_tumors.quantile([.25,.5,.75])
lowq = quartiles[0.25]
upq = quartiles[0.75]
iqr = upq-lowq
print(f"The lower quartile of Infubinol tumors is: {lowq}")
print(f"The upper quartile of Infubinol tumors is: {upq}")
print(f"The interquartile range of Infubinol tumors is: {iqr}")
lower_bound = lowq - (1.5*iqr)
upper_bound = upq + (1.5*iqr)
print(f"Numbers less than {lower_bound} might be outliers.")
print(f"Numbers more than {upper_bound} might be outliers.")


# In[44]:


#Ceftamin Data
Ceftamin_last = Ceftamin.groupby('Mouse ID').max()['Timepoint']
Ceftamin_vol = pd.DataFrame(Ceftamin_last)
Ceftamin_merge = pd.merge(Ceftamin, df, on=("Mouse ID","Timepoint"),how="left")
Ceftamin_merge.head()


# In[45]:


Ceftamin_tumors = Ceftamin_merge["Tumor Volume (mm3)_x"]
quartiles =Ceftamin_tumors.quantile([.25,.5,.75])
lowq = quartiles[0.25]
upq = quartiles[0.75]
iqr = upq-lowq
print(f"The lower quartile of Ceftamin tumors is: {lowq}")
print(f"The upper quartile of Ceftamin tumors is: {upq}")
print(f"The interquartile range of Ceftamin tumors is: {iqr}")
lower_bound = lowq - (1.5*iqr)
upper_bound = upq + (1.5*iqr)
print(f"Numbers less than {lower_bound} might be outliers.")
print(f"Numbers more than {upper_bound} might be outliers.")


# In[51]:


#Box Plots
box_data = [Capomulin_tumors, Ramicane_tumors, Infubinol_tumors, Ceftamin_tumors]
regimen = ['Capomulin', 'Ramicane', 'Infubinol','Ceftamin']
fig1, ax1 = plt.subplots()
ax1.boxplot(box_data, labels=regimen, widths = 0.8, vert=True)
ax1.set_title('Drug Regimen vs. Final Tumor Volume')
ax1.set_ylabel('Final Tumor Volume (mm3)')
ax1.set_xlabel('Drug Regimen')
plt.show()


# In[73]:


total_drugs = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]
merged_drugs = df[df["Drug Regimen"].isin(total_drugs)]
merged_drugs.head()


# In[57]:


end_time = merged_drugs.groupby(["Drug Regimen", "Mouse ID"]).agg(tumor_size=("Tumor Volume (mm3)", lambda x: x.iloc[-1]))
end_time = end_time.stack(level=0).unstack(level=0)
for drug in total_drugs:
    print(drug)


# In[68]:


treatment = 0
for drug in total_drugs:
    quartiles = end_time[drug].quantile([.25,.5,.75]).round(2)
    lowq = quartiles[0.25].round(2)
    upq = quartiles[0.75].round(2)
    iqr = round(upq-lowq,2)
    lower_bound = round(lowq - (1.5*iqr),2)
    upper_bound = round(upq + (1.5*iqr),2)
    
    if treatment == 0:
        print(f"-------------------------------------------------------")
    print(f"The lower quartile of {drug} treatments is: {lowerq}")
    print(f"The upper quartile of {drug} treatments is: {upperq}")
    print(f"The interquartile range of {drug} treatments is: {iqr}")
    print(f"Values lower than {lower_bound} might be {drug} outliers.")
    print(f"Values higher than {upper_bound} might be {drug} outliers.")
    print(f"--------------------------------------------------------")
    treatment+=1


# In[72]:


#Boxplot
box_list = []
for drug in total_drugs:
    box_list.append(list(end_time[drug].dropna()))
plt.boxplot(box_list)
plt.xlabel("Drug Regimen")
plt.ylabel("Tumor Volume")
plt.title("Drug Regimen vs. Tumor Volume")
fig = plt.figure()
plt.show() 


# In[76]:


line = Capomulin.loc[Capomulin["Mouse ID"] == "l509",:]
line.head()


# In[80]:


#Line Plot
x_axis = line["Timepoint"]
tumsiz = line["Tumor Volume (mm3)"]
fig1, ax1 = plt.subplots()
plt.plot(x_axis, tumsiz,linewidth=2, markersize=15,marker="o",color="grey", label="Fahreneit")
plt.title('Capomulin treatmeant of mouse l509')
plt.xlabel('Timepoint (Days)')
plt.ylabel('Tumor Volume (mm3)')


# In[89]:


#Scatter Plot
fig1, ax1 = plt.subplots()
average_volume =Capomulin.groupby(['Mouse ID']).mean()
plt.scatter(average_volume['Weight (g)'],average_volume['Tumor Volume (mm3)'], color="red")
plt.title('Mouse Weight vs. Average Tumor Volume')
plt.xlabel('Mouse Weight')
plt.ylabel('Averag Tumor Volume (mm3)')


# In[90]:


#Correlation and Regression
correlation = st.pearsonr(average_volume['Weight (g)'],average_volume['Tumor Volume (mm3)'])
print(f"The correlation coefficient between the mouse weight and the average tumor volume is {round(correlation[0],2)}")


# In[92]:


#Linear Regression
(slope, intercept,rvalue, pvalue, stderr)= linregress(average_volume["Weight (g)"],average_volume["Tumor Volume (mm3)"])
regression_values=average_volume["Weight (g)"]* slope + intercept
line= f"y = {round(slope, 2)} x + {round(intercept, 2)}"
plt.scatter(average_volume["Weight (g)"],average_volume["Tumor Volume (mm3)"])
plt.plot(average_volume["Weight (g)"], regression_values)
plt.annotate(line,(20,36))
plt.title("Mouse Weight vs Average Tumor Volume for Capomulin")
plt.xlabel("Mouse Weight")
plt.ylabel("Average Tumor Volume (mm3)")
print(f"R-squared is: {round(rvalue**2,3)}")
plt.show()


# In[ ]:




