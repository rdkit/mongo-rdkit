from pathlib import Path
import sys
import os
import mongomock

sys.path.append(Path('.').resolve().parent.parent)


import Database.write


def test_write():
    client = mongomock.MongoClient()
    db = client.db


