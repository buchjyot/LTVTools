# Testing Strategy and Best Practices
* Try to test most of the source code to catch potential bug or to protect feature
* Run **runAllTests** before pushing any change to source, test or demo code
* Follow **STANDARD_TEST_TEMPLATE.m** for writing a test
* Writing small unittest methods is recommended 
* Add numerical verification and reproducibility tests
* Add performance tests to moniter performance and regression in code

# runAllTests
* Runs all the available tests in this directory by going one level deeper 
for all the sub directories

# runall
* Every directory individually provides **runall** function to run all the 
tests within that directory
