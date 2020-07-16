import pytest
import pymongo
import pymongo.errors


def pytest_addoption(parser):
    parser.addoption(
        "--server", action="store", default=None, help="enter a MongoURI or 'local' to use localhost."
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "mongo: mark test as needing mongo to run")


@pytest.fixture
def mongoURI(request):
    return request.config.getoption("--server")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--server") == "local":
        try:
            pymongo.MongoClient(serverSelectionTimeoutMS=100).server_info()
            return
        except pymongo.errors.ServerSelectionTimeoutError:
            skip_mongo = pytest.mark.skip(reason="needs Mongo instance to run")
            for item in items:
                if "mongo" in item.keywords:
                    item.add_marker(skip_mongo)
    elif config.getoption("--server") and config.getoption("--server") != "local":
        try:
            pymongo.MongoClient(config.getoption("--server"), serverSelectionTimeoutMS=100).server_info()
            return
        except pymongo.errors.ServerSelectionTimeoutError:
            skip_mongo = pytest.mark.skip(reason="needs Mongo instance to run")
            for item in items:
                if "mongo" in item.keywords:
                    item.add_marker(skip_mongo)
    skip_mongo = pytest.mark.skip(reason="needs Mongo instance to run")
    for item in items:
        if "mongo" in item.keywords:
            item.add_marker(skip_mongo)

