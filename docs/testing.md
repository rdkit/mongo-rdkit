# Testing mongo-rdkit
Tests for this integration are built using the [pytest](https://docs.pytest.org/en/stable/contents.html)
testing framework. If you didn't install pytest as part of environment setup, 
you will have to do so in order to run tests:
```
pip install -U pytest
```
## Running Tests
While you are in the top-level directory (`/your path here/mongo-rdkit`),
you can run all tests simply by opening a command line and running:
```
pytest
```
Directory position is important because as of 7/10/20, 
the tests have explicit file paths to test data in the `/data` directory.

By only running `pytest`, some tests will be skipped. This is because
some tests—for instance, those that involve aggregation pipelines—only work with a functional 
MongoDB instance. In order to successfully run these tests, we have to pass in
a MongoDB URI to the testing framework:
```
pytest --server="YOUR_MONGO_URI"
```
Passing in `"local"` will run the tests on your locally hosted MongoDB
instance. 

Pytest also offers additional flexibility in running tests. For example, you can 
run a particular test file: 
```
pytest test_similarity.py
```
You can find a full list of options by running `pytest --help` or viewing the pytest documentation, 
linked above.
## Test Directory Structure
Tests for each module are located in a `test` directory at the same level as the module. 
```
mongordkit
    | Database
        | tests
            | test_create.py
            | test_write.py
    | Search 
        | tests
            | test_similarity.py
            | test_substructure.py
```

Configurations for the tests are located in `conftest.py`.



