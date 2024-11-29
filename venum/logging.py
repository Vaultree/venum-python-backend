import logging

logger = logging.getLogger('venum')
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s] %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.NOTSET)
