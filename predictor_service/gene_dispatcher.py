from threading import Thread
import time
from gene_predictor import predict_all

def process(message):
    print(message)
    def wrap():
        print(message)
        predict_all(message)

    t = Thread(target=wrap)
    t.start()
