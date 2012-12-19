Interest Rates Models
===================

This repository holds plugins implementing several interest rate models:

-Hull And White One Factor: A no-arbitrage model, which is the industry standard for modeling the future interest rate dynamic.

-Hull And White Two Factors: A no-arbitrage model, which is the industry standard for modeling the future interest rate dynamic which uses two factors in order to handle better situations like pricing a derivative whose payoff depends on rates at different maturities.

-Pelsser's squared gaussian model: A no-arbitrage model which ensures positive interest rates realizations.

-Cox-Ingersoll-Ross model (one and two factors): A a simple mean reverting short rate model with the feature to generate only positive rate. This plug-in allows you to simulate and calibrate the model on cap prices.
