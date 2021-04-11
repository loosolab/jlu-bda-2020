import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { GraphHomeComponent } from './graph-home/graph-home.component';
import { GraphTfComponent } from './graph-tf/graph-tf.component';

const routes: Routes = [
  {path: "graph_home", component: GraphHomeComponent},
  {path: "graph_tf", component: GraphTfComponent},


  {path: "", redirectTo: "graph_home", pathMatch: "full"}
];

@NgModule({
  imports: [RouterModule.forRoot(routes)],
  exports: [RouterModule]
})
export class AppRoutingModule { }
