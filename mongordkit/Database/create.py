import pymongo
from pymongo import MongoClient
from .utils import *

def createFromHostPort(dbname, host=None, port=None):
    """
    Creates a database called dbname to the desired
    mongod instance through MongoClient. By default, this
    creates a 'registration' collection and writes the standard
    registration model into it.
    :param dbname: The desired database name as a string.
    :param host: The host name.
    :param port: The port number as a string.
    :return: The database created.
    """
    client = MongoClient(host, port)
    db = client[dbname]
    collection = db['registration']
    if db.registration.count_documents({'name': 'STANDARD'}, limit=1) == 0:
        collection.insert_one(STANDARD_SETTING)
    return db

def createFromURL(dbname, url=None):
    """
    Identical to createFromHostPort, except using a url.
    :param dbname: The desired database name as a string.
    :param url: The url as a string.
    :return: The database created.
    """
    client = MongoClient(url)
    db = client[dbname]
    collection = db.registration
    if db.registration.count_documents({'name': 'STANDARD'}, limit=1) == 0:
        collection.insert_one(STANDARD_SETTING)
    return db