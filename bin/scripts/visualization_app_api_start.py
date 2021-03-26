# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 16:32:30 2021

@author: Spyro
"""

# import main Flask class and request object
from flask import Flask, request, render_template
from flask_cors import CORS, cross_origin
import api_calls

# create the Flask app
app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = ['Content-Type']
app.config['CORS_ORIGINS'] = ['*']
app.config['CORS_METHODS'] = ['GET', 'POST']

@app.route('/getTreeData', methods=["GET"])
def getTreeData():
    #result_dict_cp = {"data":[{'item': 'GM12878', 'type': '', 'belongsTo': '', 'checked': False, 'children': [{'item': 'ATF2_ENCFF210HTZ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': False, 'children': []}, {'item': 'ATF2_ENCFF210HTZ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': False, 'children': []}, {'item': 'ATF2_ENCFF210HTZ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': False, 'children': []}, {'item': 'ATF2_ENCFF210HTZ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': False, 'children': []}, {'item': 'ATF2_ENCFF210HTZ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': False, 'children': []}, {'item': 'ATF2_ENCFF210HTZ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': False, 'children': []}, {'item': 'ATF2_ENCFF210HTZ', 'type': 'tf', 'belongsTo': 'GM12878', 'checked': False, 'children': []}]}]}
    return api_calls.get_biosource_list_for_tree()
    #return result_dict_cp

# =============================================================================
# @app.route('/getGraphPaths', methods=["POST"])
# def getGraphPaths():
#     request_data = request.get_json()
#     #data = {"data": request_data, "hey": "there"}
#     #checked_dict = api_calls.getChecked(request_data)
#     
#     return api_calls.getGraphs(api_calls.getChecked(request_data))
# =============================================================================

@app.route('/getRawData', methods=["POST"])
@cross_origin()
def getRawData():
    request_data = request.get_json()
    return api_calls.getRawData(api_calls.getChecked(request_data))


if __name__ == '__main__':
    # run app in debug mode on port 5000
    app.run(debug=True, port=5000)
