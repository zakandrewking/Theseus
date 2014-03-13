import tornado.ioloop
import tornado.web
import tornado.escape
from tornado.options import define, options, parse_command_line
import json
import cPickle as pickle

from theseus import models

# define port
define("port", default=9091, type=int)

def main():
    parse_command_line()
    application.listen(options.port)
    try:
        tornado.ioloop.IOLoop.instance().start()
    except KeyboardInterrupt:
        print "bye!"

class ModelHandler(tornado.web.RequestHandler):
    def get(self, model_name):
        try:
            id_style = self.get_arguments("id_style")[0]
            model = models.load_model(model_name, id_style=id_style)
        except IndexError:
            model = models.load_model(model_name)
        # self.set_header ('Content-Type', 'text/csv')
        self.set_header('Content-Disposition', 'attachment; filename=%s.pickle' % model_name)
        self.write(pickle.dumps(model))
        self.finish()

class GetModelsHandler(tornado.web.RequestHandler):
    def get(self, path):
        data = json.dumps({'models': models.get_model_list()})
        self.write(data)
        self.finish()
                
settings = {
        "debug": "True",
        }

application = tornado.web.Application([
    (r"/models/(.*)", ModelHandler),
    (r"/models()", GetModelsHandler)
], **settings)

if __name__=="__main__":
    main()
