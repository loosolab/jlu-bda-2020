import { Component, OnInit } from '@angular/core';
import { NavigationExtras, Router } from '@angular/router';
import { FlaskApiService } from '../service/flask-api.service';

@Component({
  selector: 'app-graph-tf',
  templateUrl: './graph-tf.component.html',
  styleUrls: ['./graph-tf.component.scss']
})
export class GraphTfComponent implements OnInit {
  graph_list: any
  graph_list_names_tf: any
  graph_list_names_bio: any
  title = 'dynamic-plots';

  constructor(
    private router: Router,
    private api_service: FlaskApiService
  ) {

    console.log("starting to create graphs")
    this.createGraphs()
  }

  createGraphs() {
    //This function creates objects for graphs and tables and sorts them into arrays. The biosources and tfs are although sorted.
    var rawData: any
    rawData = this.api_service.RawGraphData.value
    let temp_bio = []
    let temp_bio_names = []

    let temp_tf = []
    let temp_tf_names_outer = []
    let temp_tf_names = []

    for (var biosource in rawData) {
      console.log(biosource)
      temp_tf = []

      for (var tf in rawData[biosource]) {
        console.log(tf)
        
        let temp_graph = []
        //create densityscatter plot object
        var temp_density_scatter = {
          data: [
            {
              x: rawData[biosource][tf]["rawData"][0][0],
              y: rawData[biosource][tf]["rawData"][0][1],
              type: 'scattergl',
              mode: "markers",
              name: "points",
              marker: {
                color: "rgb(246,19,238)",
                size: 1.5,
                opacity: 0.3
              },
            },
            {
              x: rawData[biosource][tf]["rawData"][0][0],
              y: rawData[biosource][tf]["rawData"][0][1],
              name: 'density',
              ncontours: 20,
              colorscale: 'Jet',
              reversescale: false,
              showscale: false,
              type: 'histogram2dcontour'
            },
            {
              x: rawData[biosource][tf]["rawData"][0][0],
              name: 'ATAC density',
              marker: { color: 'rgb(19,163,246)' },
              yaxis: 'y2',
              type: 'histogram'
            },
            {
              y: rawData[biosource][tf]["rawData"][0][1],
              name: 'CHIP density',
              marker: { color: 'rgb(19,163,246)' },
              xaxis: 'x2',
              type: 'histogram'
            }

          ],
          layout: {
            showlegend: false,
            autosize: false,
            width: 600,
            height: 550,
            margin: { t: 50 },
            hovermode: 'closest',
            bargap: 0,
            xaxis: {
              title: "ATAC",
              domain: [0, 0.85],
              showgrid: false,
              zeroline: false
            },
            yaxis: {
              title: "CHIP",
              domain: [0, 0.85],
              showgrid: false,
              zeroline: false
            },
            xaxis2: {
              domain: [0.85, 1],
              showgrid: false,
              zeroline: false
            },
            yaxis2: {
              domain: [0.85, 1],
              showgrid: false,
              zeroline: false
            }
          }
        };
        
        console.log(temp_density_scatter)
        //create table for weights
        var table = {
          data: [
            {
              type: "table",
              columnwidth: 20,
              header: {
                values: [["<b>ATAC</b>"], ["<b>CHIP</b>"], ["<b>Weight</b>"]],
                align: "center",
                line: { width: 1, color: 'black' },
                fill: { color: "grey" },
                font: { family: "Arial", size: 14, color: "white" }
              },
              cells: {
                values: rawData[biosource][tf]["analysisData"],
                align: "center",
                line: { color: "black", width: 1 },
                font: { family: "Arial", size: 12, color: ["black"] }

              }
            }]

        }
        //save densityscatter and table in list for each tf
        temp_graph.push(temp_density_scatter)
        temp_graph.push(table)
        //push chrs and metadata as table?
        temp_tf.push(temp_graph)
        temp_tf_names.push(tf)
      }
      
      temp_bio.push(temp_tf)
      temp_tf_names_outer.push(temp_tf_names)
      temp_tf_names = []
      temp_bio_names.push(biosource)
    }
    //set global variables on the value of local variables
    this.graph_list = temp_bio
    this.graph_list_names_tf = temp_tf_names_outer
    this.graph_list_names_bio = temp_bio_names
    console.log(this.graph_list)

    console.log(this.graph_list_names_bio)
    console.log(this.graph_list_names_tf)
    console.log(this.graph_list)
  }

  ngOnInit(): void {
  }

  goBack(){
    this.router.navigate(["/graph_home"])
  }

}
