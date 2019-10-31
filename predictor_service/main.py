from async_consumer import ReconnectingConsumer
from gene_dispatcher import process
import logging

LOG_FORMAT = ('%(levelname) -10s %(asctime)s %(name) -30s %(funcName) '
              '-35s %(lineno) -5d: %(message)s')
LOGGER = logging.getLogger(__name__)

# logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)

amqp_url = 'amqp://guest:guest@localhost:5672/%2F'

consumer = ReconnectingConsumer(amqp_url, process)
consumer.run()