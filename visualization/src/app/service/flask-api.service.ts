import { Injectable } from '@angular/core';
import { BehaviorSubject } from 'rxjs';

@Injectable({
  providedIn: 'root'
})

export class TFItemNode {
  item: string;
  type: string;
  belongsTo: string;
  checked: boolean;
  children: TFItemNode[];
}

export class FlaskApiService {
  public Tree_Data: BehaviorSubject<TFItemNode[]> = new BehaviorSubject<TFItemNode[]>([])
  public Viszalization_Data = new BehaviorSubject(Array)
  public RawGraphData = new BehaviorSubject(Array)
  private api_adress = "http://localhost:5000/"
  constructor() { }

  getTreeDataFromAPI() {
    return new Promise(resolve => {
      var rawFile = new XMLHttpRequest();
      rawFile.open("GET", this.api_adress+"getTreeData", false)
      rawFile.onload = function () {
        var res = rawFile.response
        var resolvedJSON = JSON.parse(res)
        console.log(resolvedJSON)
        resolve(resolvedJSON)
      };
      console.log("json2", rawFile)
      rawFile.send(null)
    })
  }

  setTreeData() {
    return new Promise(resolve => {
      this.getTreeDataFromAPI().then((dataobject: any) => {
        resolve(this.Tree_Data.next(dataobject.data))
      })
    })
  }

  getPathList(){
    return new Promise(resolve => {
      var xmlhttp = new XMLHttpRequest();
      xmlhttp.open("POST", this.api_adress + "getGraphPaths")
      xmlhttp.setRequestHeader("Content-type", "application/json;charset=UTF-8")
      xmlhttp.send(JSON.stringify(this.Tree_Data.value))
      xmlhttp.onload = function () {
        var res = xmlhttp.response
        var resolvedJSON = JSON.parse(res)
        console.log(resolvedJSON)
        resolve(resolvedJSON)
      };
      console.log("json2", xmlhttp)
    })
}
  setPathList(){
    return new Promise(resolve => {
      this.getPathList().then((dataobject: any) => {
        resolve(this.Viszalization_Data.next(dataobject.data))
      })
    })
  }

  getRawData(){
    return new Promise(resolve => {
      var xmlhttp = new XMLHttpRequest();
      xmlhttp.open("POST", this.api_adress +"getRawData")
      xmlhttp.setRequestHeader("Content-type", "application/json;charset=UTF-8")
      xmlhttp.send(JSON.stringify(this.Tree_Data.value))
      xmlhttp.onload = function () {
        var res = xmlhttp.response
        var resolvedJSON = JSON.parse(res)
        console.log(resolvedJSON)
        resolve(resolvedJSON)
      };
      console.log("json2", xmlhttp)
    })
  }

  setRawData(){
    return new Promise(resolve => {
      this.getRawData().then((dataobject: any) => {
        resolve(this.RawGraphData.next(dataobject.data))
      })
    })
  }
}

