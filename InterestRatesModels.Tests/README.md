Interest Rates Models Tests
===================

These are several unit tests covering all the projects in the specific repository.
Each project has its own folder containing all the sources of the tests.
The unit tests have been packaged in order to allow running them using nunit.

In order for these tests to work correctly its needed that some options are set in 
nunit, which aren't set by default:
* Test loader -> Assembly isolation
    * Run tests in a separate process per assembly *must* be selected
* Test loader -> Advanced
    * Enable Shadow Copy *must* be disabled
