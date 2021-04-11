import { Component, OnInit } from '@angular/core';
import { ActivatedRoute, Router } from '@angular/router';
import { MatButtonModule } from '@angular/material/button';
import { HomeComponent } from '../home/home.component';
import { NestedTreeControl } from '@angular/cdk/tree';
import { MatTreeFlattener, MatTreeFlatDataSource, MatTreeNestedDataSource } from '@angular/material/tree';
import { BehaviorSubject, Observable, of as observableOf } from 'rxjs';
import { FlaskApiService, TFItemNode } from '../service/flask-api.service';



const Tree_data = {
  Biosource: {
    "one": ["tf1", "tf2"],
    "two": ["tf1", "tf2"]
  }
}



@Component({
  selector: 'app-graph-home',
  templateUrl: './graph-home.component.html',
  styleUrls: ['./graph-home.component.scss'],
})
export class GraphHomeComponent implements OnInit {
  //create TreeControl classes
  nestedTreeControl: NestedTreeControl<TFItemNode>
  nestedDataSource: MatTreeNestedDataSource<TFItemNode>
  dataChange: BehaviorSubject<TFItemNode[]> = new BehaviorSubject<TFItemNode[]>([])

  filelist: any
  graphData: any

  constructor(
    private route: ActivatedRoute,
    private router: Router,
    private api_service: FlaskApiService
  ) {
    this.api_service.setTreeData().then(() => {
      console.log("from service", this.api_service.Tree_Data.value)
      console.log(this.api_service.Tree_Data.value)
      //initialize TreeControl classes
      this.nestedTreeControl = new NestedTreeControl<TFItemNode>(this._getChildren)
      this.nestedDataSource = new MatTreeNestedDataSource()
      //get TreeData from Api
      this.api_service.Tree_Data.subscribe(data => this.nestedDataSource.data = data)
      
    })
  }

  private _getChildren = (node: TFItemNode) => { return observableOf(node.children) }

  hasNestedChild = (_: number, nodeData: TFItemNode) => { return !(nodeData.type) }


  ngOnInit(): void {

  }


  goToResults() {
    //This Function sends an API Call and navigates to the next page when it got an return.
    this.api_service.setRawData().then(() =>{
      console.log("on graph home recieved data",this.api_service.RawGraphData.value)
      this.router.navigate(["/graph_tf"])
    })
    
  }

  updateAllChecked(node: TFItemNode) {
    //This Function checks if all tfs that belong to a biosource are checked
    console.log(node)
    this.nestedDataSource.data.forEach(element => {
      if (element.item == node.belongsTo) {
        element.checked = element.children.every(t => t.checked == true)
      }
    })
  }

  checkAll(node: TFItemNode) {
    //This Function selects all tfs that belong to that biosource
    console.log(node)
    if (node.checked) {
      node.children.forEach(element => {
        element.checked = true
        console.log(element)
      });
    } else {
      node.children.forEach(element => {
        element.checked = false
        console.log(element)
      });
    }
  }

  someChecked(node: TFItemNode): boolean {
    //This Function checks if only some tfs that belong to a biosource are checked
    return node.children.filter(tf => tf.checked).length > 0 && !node.checked
  }

  goBack(){
    //This Function navigates back to the start page
    this.router.navigate(["/home"])
  }
}
