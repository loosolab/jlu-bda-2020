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
  // Bar Chart
  /*graph_density = {
    data: [
      { 
        x: [1,2,3,4,1,9,9,2,5,7,8,1,9,9,9,9,9,9,9,9],
        y: [2,3,4,4,8,6,8,2,0,4,7,9,9,9,9,9,9,9,9,9],        
        type: 'scatter', 
        mode:"markers", 
        name:"points",
        marker: {
          color: "rgb(150,0,0)",
          size: 5,
          opacity: 0.5
        },  
      },
      {
        x: [1,2,3,4,1,9,9,2,5,7,8,1,9,9,9,9,9,9,9,9],
        y: [2,3,4,4,8,6,8,2,0,4,7,9,9,9,9,9,9,9,9,9],
        name: 'density',
        ncontours: 20,
        colorscale: 'Hot',
        reversescale: true,
        showscale: false,
        type: 'histogram2dcontour'
      }
    ],
    layout: {title: 'Some Data to Hover Over'}
  };
  // Line chart
  /*graph2 = {
    data: [
      { x: [1, 2, 3, 4, 5], y: [1, 4, 9, 4, 1], type: 'scatter' },
      { x: [1, 2, 3, 4, 5], y: [1, 3, 6, 9, 6], type: 'scatter' },
      { x: [1, 2, 3, 4, 5], y: [1, 2, 4, 5, 6], type: 'scatter' },
    ],
    layout: {title: 'Some Data to Highlight'}
  };*/
  /*density_scatter = {
    data: [
    { 
        x: [],
        y: [],        
        type: 'scatter', 
        mode:"markers", 
        name:"points",
        marker: {
          color: "rgb(150,0,0)",
          size: 5,
          opacity: 0.5
        },  
      },
      {
        x: [],
        y: [],
        name: 'density',
        ncontours: 20,
        colorscale: 'Hot',
        reversescale: true,
        showscale: false,
        type: 'histogram2dcontour'
      }
    ],
    layout: {title: 'Some Data to Hover Over'}
  };*/
  constructor(
    private router: Router,
    private api_service: FlaskApiService
  ) {
    //let temp_list = []
    //temp_list.push(this.graph_density)
    //temp_list.push(this.graph2)
    //this.graph_list=temp_list
    console.log("starting to create graphs")
    this.createGraphs()
  }

  createGraphs() {
    var rawData: any
    rawData = this.api_service.RawGraphData.value
    let temp_bio = []
    let temp_bio_names = []

    let temp_tf = []
    let temp_tf_names = []
    for (var biosource in rawData) {
      console.log(biosource)
      temp_tf = []
      for (var tf in rawData[biosource]) {
        console.log(tf)
        // create 3 dataobjects because of 3 plots per tf
        let temp_graph = []
        //densityscatter
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
        /*var temp_3dsurface = {
          data: [
            {
              x: rawData[biosource][tf][0][0],
              y: rawData[biosource][tf][0][1],
              z: rawData[biosource][tf][0][2],
              type: "surface"
            }
          ],
          layout: {
            title: 'Surface Plot',
          autosize: false,
          width: 500,
          height: 500,
          margin: {
            l: 65,
            r: 50,
            b: 65,
            t: 90,
          }}
        };
      var temp_contourPlot = {
        data: [
          {
            x: rawData[biosource][tf][0][0],
            y: rawData[biosource][tf][0][1],
            //z: rawData[biosource][tf][0][2],
            ncontours: 20,
            colorscale: 'Jet',
            reversescale: true,
            showscale: false,
            type: 'contour'
          }
        ],
        layout: {title: 'Contour Plot',
        xaxis: {
          range: [0,100]
        },
        yaxsis: {
          range: [0,100]
        }}
      };*/
        console.log(temp_density_scatter)
        //console.log(temp_3dsurface)
        //console.log(temp_contourPlot)
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
        console.log(table)
        //temp_graph.push(temp_density_scatter)
        //temp_graph.push(temp_contourPlot)
        //temp_graph.push(temp_3dsurface)
        temp_graph.push(temp_density_scatter)
        temp_graph.push(table)
        //temp_graph.push(temp_contourPlot)
        //temp_graph.push(temp_3dsurface)
        temp_tf.push(temp_graph)
        temp_tf_names.push(tf)

      }
      temp_bio.push(temp_tf)
      temp_bio_names.push(biosource)
    }

    this.graph_list = temp_bio
    this.graph_list_names_tf = temp_tf_names
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
