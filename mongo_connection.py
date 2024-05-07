import pymongo

# Establish a connection to the MongoDB server
client = pymongo.MongoClient('mongodb://localhost:27017/')
db = client['genomics_db']
sequences_collection = db['sequences']
