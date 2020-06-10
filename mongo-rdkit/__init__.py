#from .write import *

def createDB(client, name):
    """Given a MongoDB client CLIENT, returns a new MongoDB database with name NAME.
    It is good practice to assign this to a variable called 'db.' NAME should be a string.
    """
    return client[name]