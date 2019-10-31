from threading import Thread
import time


def dummy_slow(gene):
    time.sleep(10)
    print(gene, time.time())


def process(message):
    def wrap():
        dummy_slow(message)

    t = Thread(target=wrap)
    t.start()
