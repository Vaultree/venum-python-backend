import logging

logger = logging.getLogger('venum')  # .addHandler(logging.NullHandler())
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s ] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.DEBUG)
