# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:08:09 2026

@Jinane
"""
import pandas as pd
import pandas_datareader.data as web
import datetime

# U.S. Inflation (CPI) – Monetary Policy Regimes
# Dataset: Consumer Price Index (Monthly)
# Monthly since 1913
# 1300 observations
# Clear regime shifts (1970s inflation, 2008 crisis, COVID spike)
# Direct monetary policy implications
# Source: FRED (Federal Reserve Economic Data)
# Series ID: CPIAUCSL

start = datetime.datetime(1950, 1, 1)
end = datetime.datetime(2025, 1, 1)

cpi = web.DataReader('CPIAUCSL', 'fred', start, end)
cpi = cpi.dropna()

# convert to inflation rate
inflation = cpi.pct_change().dropna()
inflation.head()

inflation.to_csv('C:/Users/User/Documents/Records/data/test/inflation_series.csv')

# How You Can Classify
# Option A:
# High inflation regime
# Low inflation regime
# Option B:
# Pre-Volcker (pre-1980)
# Great Moderation (1985–2007)
# Post-COVID regime
# Your model can detect structural regime differences.


###############################
# CO₂ Concentration – Climate Acceleration Detection
# Dataset: Mauna Loa CO₂ (Monthly)
# Monthly since 1958 (~800 observations)
# Strong trend + seasonality
# Policy relevance: climate change mitigation
# Tests model ability to detect deterministic trend vs structural shifts
# 🔹 Classification Idea
# Linear trend vs accelerating regime
# Pre-1990 vs Post-1990 acceleration
# Seasonal vs non-seasonal components
# Very interesting for structural dynamics.

url = "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.txt"
co2 = pd.read_csv(url, delim_whitespace=True, comment='#',
                  names=["year","month","decimal","average","deseasonalized","days","std","unc"])

co2 = co2[co2["average"] > 0]
co2["date"] = pd.to_datetime(dict(year=co2.year, month=co2.month, day=1))
co2 = co2.set_index("date")

co2_series = co2["average"]
co2_series.head()

co2_series.to_csv('C:/Users/User/Documents/Records/data/test/co2_series.csv')


# U.S. Unemployment Rate – Business Cycle Classification
# Dataset: UNRATE (Monthly)
# Monthly since 1948 (~900 obs)
# Clear recession spikes
# Direct macro policy implications
# Excellent for regime detection models
# Source : FRED
# Series ID: UNRATE
# Classification Strategy
# Recession vs expansion periods
# Pre vs post crisis structural change
# High-volatility vs stable periods
# You can align with NBER recession dates.

unrate = web.DataReader('UNRATE', 'fred', start, end)
unrate = unrate.dropna()
unrate.head()

unrate.to_csv('C:/Users/User/Documents/Records/data/test/unrate_series.csv')


# 🚀 Best For Your Paper?

# If your work is theoretical (record theory + ML classification):

# Goal	 ||                 Best Dataset
# Regime detection	        Unemployment
# Trend classification	    CO₂
# Structural breaks	        Inflation
# Policy discussion depth	Inflation + Unemployment