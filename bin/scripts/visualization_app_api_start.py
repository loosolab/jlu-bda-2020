# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 16:32:30 2021

@author: Spyro
"""

# import main Flask class and request object
from flask import Flask, request, render_template
from flask_cors import CORS, cross_origin
import api_calls

# create the Flask app and configure Cors
app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = ['Content-Type']
app.config['CORS_ORIGINS'] = ['*']
app.config['CORS_METHODS'] = ['GET', 'POST']

@app.route('/getTreeData', methods=["GET"])
def getTreeData():
    #This function return a JSON object containing the structure of the tree object for the angular visualization
    return api_calls.get_biosource_list_for_tree()


@app.route('/getRawData', methods=["POST"])
@cross_origin()
def getRawData():
    #This function gets a JSON object containing the biosources and tfs the user wants to see and returns the data of these as a JSON object.
    request_data = request.get_json()
    return api_calls.getRawData(api_calls.getChecked(request_data))


if __name__ == '__main__':
    # run app in debug mode on port 5000
    app.run(debug=True, port=5000)
