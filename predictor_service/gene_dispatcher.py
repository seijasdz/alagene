from threading import Thread
import json
from gene_predictor import predict_all
from async_sender import connect

send_message = connect('amqp://guest:guest@127.0.0.1:5672/%2F?connection_attempts=3&heartbeat=3600')


def process(message):
    def wrap():
        decoded_message = json.loads(message.decode('utf-8'))
        path = predict_all(decoded_message['dna'])
        response_message = {
            '_id': decoded_message['_id'],
            'result': path
        }
        send_message(response_message)
    t = Thread(target=wrap)
    t.start()
